"""Gaussian cube file input and output utilities."""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np

from .constants import BOHR_PER_ANGSTROM
from .models import Cube


def _parse_atoms_header(line: str) -> Tuple[int, float, float, float, Optional[int]]:
    toks = line.split()
    if len(toks) < 4:
        raise ValueError(f"Bad cube atom header line: {line!r}")
    natoms_raw = int(toks[0])
    ox, oy, oz = float(toks[1]), float(toks[2]), float(toks[3])
    nvals = int(toks[4]) if len(toks) >= 5 else None
    return natoms_raw, ox, oy, oz, nvals


def _parse_int_float4(line: str) -> Tuple[int, float, float, float]:
    toks = line.split()
    if len(toks) < 4:
        raise ValueError(f"Bad cube header line: {line!r}")
    return int(toks[0]), float(toks[1]), float(toks[2]), float(toks[3])


def read_cube(path: Path, dtype: np.dtype, cube_units: str = "bohr") -> Cube:
    """Read a Gaussian cube file.

    Parameters
    ----------
    path:
        Cube file path.
    dtype:
        NumPy dtype used for the grid values.
    cube_units:
        Unit convention used for the cube header origin and axis vectors.

        ``"bohr"`` is the correct default for cube files produced by Gaussian
        ``cubegen``.  In Gaussian
        MO cube files, a negative atom count is used to indicate the presence of
        a dataset or orbital identifier section.  It must not be interpreted as
        an Angstrom coordinate flag for this workflow.

        ``"angstrom"`` is provided only for non-Gaussian cube files whose header
        coordinates are explicitly known to be Angstrom.
    """

    if cube_units not in {"bohr", "angstrom", "auto"}:
        raise ValueError(f"cube_units must be 'bohr', 'angstrom', or 'auto', got {cube_units!r}")

    text = path.read_text(errors="ignore")
    lines = text.splitlines()
    if len(lines) < 6:
        raise ValueError(f"{path} does not look like a cube file.")

    comment1 = lines[0].rstrip("\n")
    comment2 = lines[1].rstrip("\n")

    natoms_raw, ox, oy, oz, nvals = _parse_atoms_header(lines[2])
    natoms = abs(natoms_raw)

    nx, vxx, vxy, vxz = _parse_int_float4(lines[3])
    ny, vyx, vyy, vyz = _parse_int_float4(lines[4])
    nz, vzx, vzy, vzz = _parse_int_float4(lines[5])

    # Important correction:
    # The previous implementation used natoms_raw < 0 to infer Angstrom units.
    # This is incorrect for Gaussian MO cube files.  For Gaussian cubegen output,
    # a negative atom count indicates that a dataset/orbital identifier line is
    # present after the atom block.  The origin and grid vectors remain in the
    # Gaussian cube coordinate system used by the present independent Python
    # implementation, namely Bohr for numerical moment-density integration.
    #
    # Therefore, "auto" intentionally resolves to "bohr" for this program.
    if cube_units == "angstrom":
        unit_scale_to_bohr = BOHR_PER_ANGSTROM
    else:
        unit_scale_to_bohr = 1.0

    origin_bohr = unit_scale_to_bohr * np.array([ox, oy, oz], dtype=float)
    vx_bohr = unit_scale_to_bohr * np.array([vxx, vxy, vxz], dtype=float)
    vy_bohr = unit_scale_to_bohr * np.array([vyx, vyy, vyz], dtype=float)
    vz_bohr = unit_scale_to_bohr * np.array([vzx, vzy, vzz], dtype=float)

    atom_lines = lines[6: 6 + natoms]
    if len(atom_lines) != natoms:
        raise ValueError(f"{path}: expected {natoms} atom lines, got {len(atom_lines)}")

    cursor = 6 + natoms
    val_ids: Optional[List[int]] = None
    if nvals is not None:
        if cursor >= len(lines):
            raise ValueError(f"{path}: expected dataset-id line after atoms (nvals={nvals}), but file ended.")
        id_toks = lines[cursor].split()
        if not id_toks:
            raise ValueError(f"{path}: empty dataset-id line after atoms.")
        val_ids = [int(x) for x in id_toks]
        cursor += 1

    data_text = "\n".join(lines[cursor:]).strip()
    arr = np.fromstring(data_text, sep=" ", dtype=dtype)
    expected = nx * ny * nz * (nvals if nvals is not None else 1)
    if arr.size < expected:
        raise ValueError(f"{path}: cube data short. expected {expected}, got {arr.size}.")
    if arr.size > expected:
        arr = arr[:expected]

    if nvals is None or nvals == 1:
        data = arr.reshape((nx, ny, nz), order="C").astype(dtype, copy=False)
    else:
        data = arr[: nx * ny * nz].reshape((nx, ny, nz), order="C").astype(dtype, copy=False)

    return Cube(
        comment1=comment1,
        comment2=comment2,
        natoms_raw=natoms_raw,
        natoms=natoms,
        nvals=nvals,
        val_ids=val_ids,
        origin_bohr=origin_bohr,
        nx=nx,
        ny=ny,
        nz=nz,
        vx_bohr=vx_bohr,
        vy_bohr=vy_bohr,
        vz_bohr=vz_bohr,
        atoms_block=atom_lines,
        data=data,
        unit_scale_to_bohr=unit_scale_to_bohr,
    )


def write_cube(path: Path, cube: Cube, data: np.ndarray, comment1: str, comment2: str) -> None:
    """Write a cube file while preserving the unit convention of the reference cube."""

    if data.shape != (cube.nx, cube.ny, cube.nz):
        raise ValueError(f"Data shape {data.shape} does not match cube grid {(cube.nx, cube.ny, cube.nz)}")

    scale_from_bohr = 1.0 / cube.unit_scale_to_bohr
    origin = scale_from_bohr * cube.origin_bohr
    vx = scale_from_bohr * cube.vx_bohr
    vy = scale_from_bohr * cube.vy_bohr
    vz = scale_from_bohr * cube.vz_bohr

    def to_ascii(s: str) -> str:
        return s.encode("ascii", "replace").decode("ascii")

    c1 = to_ascii(comment1)
    c2 = to_ascii(comment2)

    with path.open("w", encoding="utf-8") as f:
        f.write(f"{c1}\n")
        f.write(f"{c2}\n")
        if cube.nvals is None:
            f.write(f"{cube.natoms_raw:5d} {origin[0]: .6f} {origin[1]: .6f} {origin[2]: .6f}\n")
        else:
            f.write(
                f"{cube.natoms_raw:5d} {origin[0]: .6f} {origin[1]: .6f} {origin[2]: .6f} {cube.nvals:d}\n"
            )
        f.write(f"{cube.nx:5d} {vx[0]: .6f} {vx[1]: .6f} {vx[2]: .6f}\n")
        f.write(f"{cube.ny:5d} {vy[0]: .6f} {vy[1]: .6f} {vy[2]: .6f}\n")
        f.write(f"{cube.nz:5d} {vz[0]: .6f} {vz[1]: .6f} {vz[2]: .6f}\n")
        for ln in cube.atoms_block:
            f.write(ln.rstrip("\n") + "\n")

        if cube.nvals is not None:
            if cube.val_ids:
                f.write(" ".join(str(x) for x in cube.val_ids) + "\n")
            else:
                f.write("1 0\n")

        flat = data.reshape(-1, order="C")
        for i in range(0, flat.size, 6):
            chunk = flat[i: i + 6]
            f.write(" ".join(f"{float(v): .5E}" for v in chunk) + "\n")


def coordinate_1d_aligned_bohr(cube: Cube) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    tol = 1e-8
    vx, vy, vz = cube.vx_bohr, cube.vy_bohr, cube.vz_bohr
    aligned = (
        abs(vx[1]) < tol
        and abs(vx[2]) < tol
        and abs(vy[0]) < tol
        and abs(vy[2]) < tol
        and abs(vz[0]) < tol
        and abs(vz[1]) < tol
    )
    if not aligned:
        raise ValueError("Cube axes not aligned.")
    x = cube.origin_bohr[0] + np.arange(cube.nx, dtype=float) * vx[0]
    y = cube.origin_bohr[1] + np.arange(cube.ny, dtype=float) * vy[1]
    z = cube.origin_bohr[2] + np.arange(cube.nz, dtype=float) * vz[2]
    return x, y, z


def coordinate_3d_general_bohr(cube: Cube) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    ix = np.arange(cube.nx, dtype=float)[:, None, None]
    iy = np.arange(cube.ny, dtype=float)[None, :, None]
    iz = np.arange(cube.nz, dtype=float)[None, None, :]
    x = cube.origin_bohr[0] + ix * cube.vx_bohr[0] + iy * cube.vy_bohr[0] + iz * cube.vz_bohr[0]
    y = cube.origin_bohr[1] + ix * cube.vx_bohr[1] + iy * cube.vy_bohr[1] + iz * cube.vz_bohr[1]
    z = cube.origin_bohr[2] + ix * cube.vx_bohr[2] + iy * cube.vy_bohr[2] + iz * cube.vz_bohr[2]
    return x, y, z


def cube_step_sizes_bohr(cube: Cube) -> Tuple[float, float, float]:
    dx = float(np.linalg.norm(cube.vx_bohr))
    dy = float(np.linalg.norm(cube.vy_bohr))
    dz = float(np.linalg.norm(cube.vz_bohr))
    if dx == 0.0 or dy == 0.0 or dz == 0.0:
        raise ValueError("Invalid cube axis vectors with zero norm.")
    return dx, dy, dz


def cube_voxel_volume_bohr3(cube: Cube) -> float:
    """Return the absolute voxel volume from the three cube lattice vectors."""

    mat = np.vstack([cube.vx_bohr, cube.vy_bohr, cube.vz_bohr])
    return float(abs(np.linalg.det(mat)))
