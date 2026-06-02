"""Numerical construction of EDTM and MDTM density grids in CPLkit.

This module provides two accumulation paths.

1. Streaming accumulation, which is memory-stable but loops over TD-DFT
   configurations in Python.
2. NumPy block accumulation, which rewrites the dominant excitation loop as
   dense matrix products over blocks of grid points. This path is much faster
   for large IOp(9/40=4/5) outputs when sufficient RAM or scratch disk is
   available.
"""

from __future__ import annotations

import math
import time
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

from .cubeio import coordinate_1d_aligned_bohr, coordinate_3d_general_bohr, cube_step_sizes_bohr
from .models import Cube, Excitation
from .utils import Timer, progress_line

ProgressCallback = Optional[Callable[[str], None]]


def orbital_gradients(phi: np.ndarray, cube: Cube) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Evaluate orbital gradients by finite differences on the cube grid."""

    dx, dy, dz = cube_step_sizes_bohr(cube)
    dphidx, dphidy, dphidz = np.gradient(phi, dx, dy, dz, edge_order=2)
    return dphidx.astype(phi.dtype, copy=False), dphidy.astype(phi.dtype, copy=False), dphidz.astype(phi.dtype, copy=False)


def excitation_weight(
    coeff: float,
    density_mode: str = "contribution-map",
    *,
    validation_edtm_factor: float = -2.0,
    validation_mdtm_factor: float = -1.0,
) -> float:
    """Return the EDTM accumulation weight for an excitation coefficient.

    Two modes are intentionally separated.

    ``contribution-map`` uses a signed square visualization convention,
    ``w = 2*c*abs(c)``.  The magnitude follows the closed-shell configuration
    contribution ``2*c**2`` that is commonly read from Gaussian TD-DFT output,
    while the sign follows the printed CI coefficient.  This convention is
    useful for spatial contribution maps, but it is not a coefficient-linear
    transition-density amplitude and is not expected to reproduce Gaussian
    EDTM or MDTM components by direct integration.

    ``validation`` uses a coefficient-linear transition-density amplitude for
    comparison with Gaussian transition moments.  The default EDTM factor is
    ``-2`` so that the electron charge sign is included for the electric
    transition dipole moment density.
    """

    if density_mode == "contribution-map":
        return 2.0 * coeff * abs(coeff)
    if density_mode == "validation":
        return validation_edtm_factor * coeff
    raise ValueError(f"Unknown density_mode: {density_mode!r}")


def excitation_weights(
    coeff: float,
    density_mode: str = "contribution-map",
    *,
    validation_edtm_factor: float = -2.0,
    validation_mdtm_factor: float = -1.0,
) -> Tuple[float, float]:
    """Return ``(EDTM weight, MDTM scale)`` for one excitation.

    In ``contribution-map`` mode, this uses the signed square visualization convention
    used for contribution analysis of dominant orbital transitions:

        EDTM weight = 2*c*abs(c)
        MDTM scale  = 0.5*EDTM weight

    This weighting maps signed configuration contributions. It is deliberately
    separated from the coefficient-linear validation mode.

    In ``validation`` mode, the EDTM and MDTM prefactors are exposed as
    parameters because Gaussian's printed TD-DFT coefficients and magnetic
    transition-moment conventions may require diagnostic sign or prefactor
    checks.  The defaults implement a closed-shell singlet, coefficient-linear
    starting point:

        EDTM weight = -2*c
        MDTM scale  = -1.0*c
    """

    if density_mode == "contribution-map":
        w = 2.0 * coeff * abs(coeff)
        return w, 0.5 * w
    if density_mode == "validation":
        return validation_edtm_factor * coeff, validation_mdtm_factor * coeff
    raise ValueError(f"Unknown density_mode: {density_mode!r}")


def _aligned_coords(cube: Cube, dtype: np.dtype) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    x, y, z = coordinate_1d_aligned_bohr(cube)
    return x.astype(dtype, copy=False), y.astype(dtype, copy=False), z.astype(dtype, copy=False)


def _general_coords(cube: Cube, dtype: np.dtype) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    x, y, z = coordinate_3d_general_bohr(cube)
    return x.astype(dtype, copy=False), y.astype(dtype, copy=False), z.astype(dtype, copy=False)


def precompute_r_cross_grad(
    cube: Cube,
    grad: Tuple[np.ndarray, np.ndarray, np.ndarray],
    coords_mode: str,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute r cross grad(phi) with dtype-preserving temporary arrays.

    This function is retained for compatibility. The command-line interface uses
    either the streaming functions or the NumPy block accumulator below.
    """

    dphidx, dphidy, dphidz = grad
    dtype = dphidx.dtype
    ax = np.empty_like(dphidx, dtype=dtype)
    ay = np.empty_like(dphidx, dtype=dtype)
    az = np.empty_like(dphidx, dtype=dtype)
    tmp = np.empty_like(dphidx, dtype=dtype)

    if coords_mode == "aligned":
        x, y, z = _aligned_coords(cube, dtype)
        np.multiply(y[None, :, None], dphidz, out=ax)
        np.multiply(z[None, None, :], dphidy, out=tmp)
        ax -= tmp

        np.multiply(z[None, None, :], dphidx, out=ay)
        np.multiply(x[:, None, None], dphidz, out=tmp)
        ay -= tmp

        np.multiply(x[:, None, None], dphidy, out=az)
        np.multiply(y[None, :, None], dphidx, out=tmp)
        az -= tmp
        return ax, ay, az

    x, y, z = _general_coords(cube, dtype)
    np.multiply(y, dphidz, out=ax)
    np.multiply(z, dphidy, out=tmp)
    ax -= tmp

    np.multiply(z, dphidx, out=ay)
    np.multiply(x, dphidz, out=tmp)
    ay -= tmp

    np.multiply(x, dphidy, out=az)
    np.multiply(y, dphidx, out=tmp)
    az -= tmp
    return ax, ay, az


# -----------------------------------------------------------------------------
# Low-memory streaming path
# -----------------------------------------------------------------------------


def accumulate_edtm_components(
    cube: Cube,
    coords_mode: str,
    phi_occ: np.ndarray,
    phi_virt: np.ndarray,
    weight: float,
    rho: Tuple[np.ndarray, np.ndarray, np.ndarray],
) -> None:
    """Accumulate EDTM component densities for one excitation in a memory-stable way."""

    dtype = phi_occ.dtype
    tau = np.empty_like(phi_occ, dtype=dtype)
    tmp = np.empty_like(phi_occ, dtype=dtype)
    np.multiply(phi_occ, phi_virt, out=tau)

    if coords_mode == "aligned":
        x, y, z = _aligned_coords(cube, dtype)
        coords = (x[:, None, None], y[None, :, None], z[None, None, :])
    else:
        coords = _general_coords(cube, dtype)

    for comp, coord in enumerate(coords):
        np.multiply(coord, tau, out=tmp)
        tmp *= weight
        np.add(rho[comp], tmp, out=rho[comp])


def _a_component_into(
    out: np.ndarray,
    tmp: np.ndarray,
    coords_mode: str,
    coords: Tuple[np.ndarray, np.ndarray, np.ndarray],
    grad: Tuple[np.ndarray, np.ndarray, np.ndarray],
    component: int,
) -> None:
    """Write one component of r cross grad(phi) into out."""

    dphidx, dphidy, dphidz = grad
    x, y, z = coords

    if component == 0:
        if coords_mode == "aligned":
            np.multiply(y[None, :, None], dphidz, out=out)
            np.multiply(z[None, None, :], dphidy, out=tmp)
        else:
            np.multiply(y, dphidz, out=out)
            np.multiply(z, dphidy, out=tmp)
        out -= tmp
        return

    if component == 1:
        if coords_mode == "aligned":
            np.multiply(z[None, None, :], dphidx, out=out)
            np.multiply(x[:, None, None], dphidz, out=tmp)
        else:
            np.multiply(z, dphidx, out=out)
            np.multiply(x, dphidz, out=tmp)
        out -= tmp
        return

    if coords_mode == "aligned":
        np.multiply(x[:, None, None], dphidy, out=out)
        np.multiply(y[None, :, None], dphidx, out=tmp)
    else:
        np.multiply(x, dphidy, out=out)
        np.multiply(y, dphidx, out=tmp)
    out -= tmp


def accumulate_mdtm_components(
    cube: Cube,
    coords_mode: str,
    phi_occ: np.ndarray,
    phi_virt: np.ndarray,
    grad_occ: Tuple[np.ndarray, np.ndarray, np.ndarray],
    grad_virt: Tuple[np.ndarray, np.ndarray, np.ndarray],
    scale: float,
    rho: Tuple[np.ndarray, np.ndarray, np.ndarray],
) -> None:
    """Accumulate MDTM component densities for one excitation without caching A fields.

    The implemented expression is
        rho_m += scale * [phi_virt * (r x grad phi_occ) - phi_occ * (r x grad phi_virt)].
    Components are accumulated one at a time so that six full A-component arrays
    are not kept in memory.
    """

    dtype = phi_occ.dtype
    a = np.empty_like(phi_occ, dtype=dtype)
    tmp = np.empty_like(phi_occ, dtype=dtype)

    if coords_mode == "aligned":
        coords = _aligned_coords(cube, dtype)
    else:
        coords = _general_coords(cube, dtype)

    for comp in range(3):
        _a_component_into(a, tmp, coords_mode, coords, grad_occ, comp)
        np.multiply(a, phi_virt, out=tmp)
        tmp *= scale
        np.add(rho[comp], tmp, out=rho[comp])

        _a_component_into(a, tmp, coords_mode, coords, grad_virt, comp)
        np.multiply(a, phi_occ, out=tmp)
        tmp *= scale
        np.subtract(rho[comp], tmp, out=rho[comp])


# -----------------------------------------------------------------------------
# NumPy block path
# -----------------------------------------------------------------------------


def _weight_matrices(
    excitations: Sequence[Excitation],
    occ_mos: Sequence[int],
    virt_mos: Sequence[int],
    dtype: np.dtype,
    density_mode: str = "contribution-map",
    *,
    validation_edtm_factor: float = -2.0,
    validation_mdtm_factor: float = -1.0,
) -> Tuple[np.ndarray, np.ndarray]:
    """Build dense occ x virt coefficient matrices for EDTM and MDTM accumulation."""

    occ_index = {mo: i for i, mo in enumerate(occ_mos)}
    virt_index = {mo: i for i, mo in enumerate(virt_mos)}
    w_edtm = np.zeros((len(occ_mos), len(virt_mos)), dtype=dtype)
    w_mdtm = np.zeros_like(w_edtm)
    for ex in excitations:
        oi = occ_index[ex.occ]
        vi = virt_index[ex.virt]
        w_edtm_ex, w_mdtm_ex = excitation_weights(
            ex.coeff,
            density_mode,
            validation_edtm_factor=validation_edtm_factor,
            validation_mdtm_factor=validation_mdtm_factor,
        )
        w_edtm[oi, vi] += dtype.type(w_edtm_ex) if hasattr(dtype, "type") else w_edtm_ex
        w_mdtm[oi, vi] += dtype.type(w_mdtm_ex) if hasattr(dtype, "type") else w_mdtm_ex
    return w_edtm, w_mdtm


def _estimate_stack_gb(cube: Cube, n_mos: int, dtype: np.dtype) -> float:
    itemsize = np.dtype(dtype).itemsize
    return cube.nx * cube.ny * cube.nz * n_mos * itemsize / (1024.0**3)


def _safe_mmap_name(baseprefix: str, state_tag: str, dtype: np.dtype, n_mos: int) -> str:
    clean = "".join(ch if ch.isalnum() or ch in "._" else "_" for ch in baseprefix)
    clean = clean[:60]
    return f"{clean}_{state_tag}_mo_stack_{np.dtype(dtype).name}_{n_mos}.dat"


def build_mo_stack(
    cube: Cube,
    mos: Sequence[int],
    get_phi: Callable[[int], np.ndarray],
    dtype: np.dtype,
    stack_mode: str,
    scratch_dir: Path,
    stack_label: str,
    progress: ProgressCallback = None,
) -> np.ndarray:
    """Build an MO stack with shape (nx, ny, nz, n_mos).

    In memory mode, a NumPy array is used. In memmap mode, the stack is backed
    by a scratch file. The latter avoids allocating a very large contiguous array
    in RAM, while still enabling vectorized slab operations.
    """

    shape = (cube.nx, cube.ny, cube.nz, len(mos))
    if stack_mode == "memmap":
        scratch_dir.mkdir(parents=True, exist_ok=True)
        stack_path = scratch_dir / stack_label
        stack = np.memmap(stack_path, dtype=dtype, mode="w+", shape=shape, order="C")
        if progress:
            progress(f"[Vectorized] MO stack mode=memmap path={stack_path}")
    else:
        stack = np.empty(shape, dtype=dtype, order="C")
        if progress:
            progress("[Vectorized] MO stack mode=memory")

    t0 = time.perf_counter()
    n = len(mos)
    for j, mo in enumerate(mos, 1):
        phi = get_phi(mo).astype(dtype, copy=False)
        stack[:, :, :, j - 1] = phi
        if progress and (j == 1 or j % 10 == 0 or j == n):
            progress(progress_line("  stack", j, n, t0))
    if isinstance(stack, np.memmap):
        stack.flush()
    return stack


def _matmul_contract(phi_o: np.ndarray, weights: np.ndarray, phi_v: np.ndarray) -> np.ndarray:
    """Return row-wise phi_o @ weights dotted with phi_v."""

    tmp = phi_o @ weights
    return np.einsum("ij,ij->i", tmp, phi_v, optimize=True)


def _a_component_from_grad(
    dphidx: np.ndarray,
    dphidy: np.ndarray,
    dphidz: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    component: int,
) -> np.ndarray:
    """Return one component of r cross grad(phi) for a slab.

    The input arrays have shape (sx, ny, nz, n_mos). x, y, and z are 1D
    coordinate vectors for the same slab and the full y and z axes.
    """

    if component == 0:
        return y[None, :, None, None] * dphidz - z[None, None, :, None] * dphidy
    if component == 1:
        return z[None, None, :, None] * dphidx - x[:, None, None, None] * dphidz
    return x[:, None, None, None] * dphidy - y[None, :, None, None] * dphidx


def accumulate_vectorized_blocks(
    cube: Cube,
    coords_mode: str,
    excitations: Sequence[Excitation],
    get_phi: Callable[[int], np.ndarray],
    dtype: np.dtype,
    rho_edtm: Tuple[np.ndarray, np.ndarray, np.ndarray],
    rho_mdtm: Tuple[np.ndarray, np.ndarray, np.ndarray],
    *,
    x_block: int = 1,
    stack_mode: str = "memory",
    scratch_dir: Optional[Path] = None,
    stack_label: str = "mo_stack.dat",
    progress: ProgressCallback = None,
    density_mode: str = "contribution-map",
    validation_edtm_factor: float = -2.0,
    validation_mdtm_factor: float = -1.0,
) -> None:
    """Accumulate EDTM and MDTM by NumPy matrix products over x-slabs.

    The streaming implementation loops over excitations in Python. This function
    instead constructs coefficient matrices W and S and evaluates, for each block
    of grid points,

        tau(r) = Phi_occ(r) W Phi_virt(r)^T

    by BLAS-backed NumPy matrix multiplication. MDTM is evaluated in the same
    block using A(r) = r x grad(phi)(r):

        rho_m = A_occ S Phi_virt^T - Phi_occ S A_virt^T.

    The implementation currently requires aligned cube axes. Non-aligned grids
    should use the streaming path.
    """

    if coords_mode != "aligned":
        raise ValueError("Vectorized accumulation currently requires aligned cube axes. Use --accumulation stream.")

    if x_block < 1:
        raise ValueError("x_block must be >= 1.")

    occ_mos = sorted({ex.occ for ex in excitations})
    virt_mos = sorted({ex.virt for ex in excitations})
    all_mos = sorted(set(occ_mos) | set(virt_mos))
    all_index = {mo: i for i, mo in enumerate(all_mos)}
    occ_cols = np.array([all_index[mo] for mo in occ_mos], dtype=np.int64)
    virt_cols = np.array([all_index[mo] for mo in virt_mos], dtype=np.int64)
    w_edtm, w_mdtm = _weight_matrices(
        excitations,
        occ_mos,
        virt_mos,
        dtype,
        density_mode,
        validation_edtm_factor=validation_edtm_factor,
        validation_mdtm_factor=validation_mdtm_factor,
    )

    if progress:
        progress(
            f"[Vectorized] unique occupied={len(occ_mos)} unique virtual={len(virt_mos)} "
            f"unique total={len(all_mos)} excitations={len(excitations)}"
        )
        progress(f"[Vectorized] estimated MO stack size={_estimate_stack_gb(cube, len(all_mos), dtype):.2f} GB")

    if scratch_dir is None:
        scratch_dir = Path(".")
    stack = build_mo_stack(cube, all_mos, get_phi, dtype, stack_mode, scratch_dir, stack_label, progress)

    dx, dy, dz = cube_step_sizes_bohr(cube)
    x, y, z = _aligned_coords(cube, dtype)

    nx, ny, nz = cube.nx, cube.ny, cube.nz
    n_blocks = math.ceil(nx / x_block)
    t0 = time.perf_counter()

    # Cast weights explicitly for BLAS consistency.
    w_edtm = np.asarray(w_edtm, dtype=dtype, order="C")
    w_mdtm = np.asarray(w_mdtm, dtype=dtype, order="C")

    for iblock, x0 in enumerate(range(0, nx, x_block), 1):
        x1 = min(nx, x0 + x_block)

        # Use a padded slab so finite differences at the block edges have the
        # same neighboring grid points as the full-array gradient. The minimum
        # length of 3 is needed for edge_order=2.
        pad0 = max(0, x0 - 2)
        pad1 = min(nx, x1 + 2)
        if pad1 - pad0 < 3:
            if pad0 == 0:
                pad1 = min(nx, 3)
            else:
                pad0 = max(0, nx - 3)
        slab = np.asarray(stack[pad0:pad1, :, :, :], dtype=dtype)
        dphidx, dphidy, dphidz = np.gradient(slab, dx, dy, dz, axis=(0, 1, 2), edge_order=2)

        trim0 = x0 - pad0
        trim1 = trim0 + (x1 - x0)
        phi = slab[trim0:trim1, :, :, :]
        gx = dphidx[trim0:trim1, :, :, :].astype(dtype, copy=False)
        gy = dphidy[trim0:trim1, :, :, :].astype(dtype, copy=False)
        gz = dphidz[trim0:trim1, :, :, :].astype(dtype, copy=False)

        b = (x1 - x0) * ny * nz
        phi2 = phi.reshape((b, len(all_mos)), order="C")
        phi_o = np.asarray(phi2[:, occ_cols], dtype=dtype, order="C")
        phi_v = np.asarray(phi2[:, virt_cols], dtype=dtype, order="C")

        tau = _matmul_contract(phi_o, w_edtm, phi_v).reshape((x1 - x0, ny, nz))

        np.add(rho_edtm[0][x0:x1, :, :], x[x0:x1, None, None] * tau, out=rho_edtm[0][x0:x1, :, :])
        np.add(rho_edtm[1][x0:x1, :, :], y[None, :, None] * tau, out=rho_edtm[1][x0:x1, :, :])
        np.add(rho_edtm[2][x0:x1, :, :], z[None, None, :] * tau, out=rho_edtm[2][x0:x1, :, :])

        phi_o_s = phi_o @ w_mdtm

        x_slab = x[x0:x1]
        for comp in range(3):
            a = _a_component_from_grad(gx, gy, gz, x_slab, y, z, comp)
            a2 = a.reshape((b, len(all_mos)), order="C")
            a_o = np.asarray(a2[:, occ_cols], dtype=dtype, order="C")
            a_v = np.asarray(a2[:, virt_cols], dtype=dtype, order="C")

            term1 = np.einsum("ij,ij->i", a_o @ w_mdtm, phi_v, optimize=True)
            term2 = np.einsum("ij,ij->i", phi_o_s, a_v, optimize=True)
            delta = (term1 - term2).reshape((x1 - x0, ny, nz))
            np.add(rho_mdtm[comp][x0:x1, :, :], delta, out=rho_mdtm[comp][x0:x1, :, :])

        if progress and (iblock == 1 or iblock % 2 == 0 or iblock == n_blocks):
            progress(progress_line("  blocks", iblock, n_blocks, t0))

    if isinstance(stack, np.memmap):
        stack.flush()
