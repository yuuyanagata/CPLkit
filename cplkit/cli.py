"""Command-line interface for the independent CPLkit Python package."""

from __future__ import annotations

import argparse
import math
import time
from pathlib import Path
from collections import OrderedDict
from typing import Dict, List, Optional, Tuple

import numpy as np

from .cubeio import (
    coordinate_1d_aligned_bohr,
    coordinate_3d_general_bohr,
    read_cube,
    write_cube,
)
from .cubegen import ensure_mo_cubes
from .densities import (
    accumulate_edtm_components,
    accumulate_mdtm_components,
    accumulate_vectorized_blocks,
    excitation_weights,
    orbital_gradients,
)
from .gaussian import (
    build_cpl_rows,
    diagnose_coefficients,
    parse_excited_state_block,
    parse_edtm_vector,
    parse_mdtm_vector,
    write_cpl_csv,
)
from .models import CPLRow
from .utils import Timer, eprint, progress_line
from .validation import compare_component_integrals, write_validation_csv
from .cubeio import cube_voxel_volume_bohr3
from . import IMPLEMENTATION_NOTE


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description=(
            "Generate EDTM and MDTM "
            "density cube files from Gaussian TD-DFT output. CPLkit is an independent "
            "Python implementation and is not a MATLAB port or wrapper."
        )
    )
    ap.add_argument("--log", required=True, type=Path, help="Gaussian TD-DFT output (*.log/*.out)")
    ap.add_argument("--state", type=int, default=None, help="Excited state index (1-based) for cube generation")
    ap.add_argument("--chk", type=Path, default=None, help="Checkpoint or formatted checkpoint file for cubegen")
    ap.add_argument("--mocubes_dir", type=Path, default=None, help="Directory with mo<MO>.cube files")
    ap.add_argument("--cubegen", default="cubegen", help="Path to cubegen executable")
    ap.add_argument("--cubegen_npts", type=int, default=0, help="First cubegen argument")
    ap.add_argument("--cubegen_grid", default="-3 h", help='Remaining cubegen arguments, for example "-3 h"')
    ap.add_argument("--cubegen_prefix", default="", help="Extra options inserted after the cubegen executable and before the nprocs argument, for example: --backend cuda --threads 32 --grid fast")
    ap.add_argument("--cubegen_jobs", type=int, default=1, help="Number of concurrent cubegen processes used for MO cube generation. Use 2 to 4 on Windows unless sufficient memory and disk I/O are available.")
    ap.add_argument("--cubegen_multi_mo_chunk", type=int, default=1, help="Number of orbitals requested in each Gaussian-compatible cube generatorCpp multi-MO call. Values greater than 1 enable one call to request multiple MOs; the current Gaussian-compatible cube generatorCpp writes one scalar cube per requested orbital.")
    ap.add_argument("--keep_multi_mo_cubes", action="store_true", help="Reserved for compatibility with older Gaussian-compatible cube generatorCpp builds that produced multi-valued cubes. The current separate-output mode does not create a multi-valued cube to keep.")
    ap.add_argument("--outdir", type=Path, default=Path("."), help="Output directory")
    ap.add_argument("--outprefix", default=None, help="Output prefix. Default is S<state>")
    ap.add_argument("--overwrite_mo_cubes", action="store_true", help="Force regeneration of MO cubes")
    ap.add_argument("--keep_components", action="store_true", help="Also write x, y, z component cubes")
    ap.add_argument("--dtype", choices=["float32", "float64"], default="float32", help="Internal dtype")
    ap.add_argument("--max_phi_cache", type=int, default=4, help="Maximum number of MO cube arrays retained in memory in the streaming path. The vectorized path keeps its own MO stack.")
    ap.add_argument("--cache_a", action="store_true", help="Cache r cross grad(phi) fields. Retained for compatibility; disabled by default.")
    ap.add_argument("--coords", choices=["auto", "aligned", "general"], default="auto", help="Coordinate handling mode")
    ap.add_argument("--cube_units", choices=["bohr", "angstrom", "auto"], default="bohr", help="Unit convention for cube header origin and axis vectors. Use bohr for Gaussian cubegen output. auto is treated as bohr for Gaussian-compatible cubes.")
    ap.add_argument("--accumulation", choices=["stream", "vectorized"], default="stream", help="Use the streaming excitation-by-excitation loop or the NumPy block matrix accumulator for EDTM/MDTM reconstruction.")
    ap.add_argument("--density_mode", choices=["contribution-map", "validation"], default="contribution-map", help="Density reconstruction convention. contribution-map uses signed square 2*c*abs(c) weights to display configuration contributions. validation uses coefficient-linear amplitudes for Gaussian EDTM/MDTM transition-moment checks.")
    ap.add_argument("--phase_align", choices=["none", "edtm"], default="none", help="Optionally multiply all reconstructed densities by -1 when the integrated EDTM is antiparallel to the Gaussian EDTM. This removes arbitrary excited-state phase differences for validation tables.")
    ap.add_argument("--validation_edtm_factor", type=float, default=-2.0, help="EDTM prefactor used only in --density_mode validation. EDTM weight is validation_edtm_factor*c. Default includes the electron charge sign: -2*c.")
    ap.add_argument("--validation_mdtm_factor", type=float, default=-1.0, help="MDTM prefactor used only in --density_mode validation. MDTM scale is validation_mdtm_factor*c for [phi_v A(phi_i)-phi_i A(phi_v)]. Default is -1.0*c.")
    ap.add_argument("--vectorized_x_block", type=int, default=1, help="Number of x-grid slices processed per NumPy block. Use 1 or 2 for large IOp(9/40=5) jobs.")
    ap.add_argument("--vectorized_stack_mode", choices=["memory", "memmap"], default="memmap", help="Storage for the MO stack used by --accumulation vectorized. memmap is safer on Windows; memory is faster if RAM is sufficient.")
    ap.add_argument("--vectorized_scratch_dir", type=Path, default=None, help="Scratch directory for the vectorized memmap MO stack. Default is <outdir>/_scratch.")
    ap.add_argument("--require_vectorized", action="store_true", help="Abort if --accumulation vectorized cannot be used. This prevents silent fallback to the slower streaming path.")
    ap.add_argument("--cpl_only", action="store_true", help="Only export the CPL CSV and skip cube generation")
    ap.add_argument("--no_cpl_csv", action="store_true", help="Disable CPL CSV export")
    ap.add_argument("--cpl_csv_path", type=Path, default=None, help="Output path for the CPL CSV")
    ap.add_argument(
        "--strict_ci",
        action="store_true",
        help="Abort cube generation when the printed CI coefficients appear truncated.",
    )
    ap.add_argument(
        "--skip_validation_report",
        action="store_true",
        help="Do not write cube-integral validation diagnostics for density reconstruction.",
    )
    ap.add_argument(
        "--validation_csv_path",
        type=Path,
        default=None,
        help="Output path for the moment-density validation CSV.",
    )
    ap.add_argument(
        "--fail_on_validation",
        action="store_true",
        help="Abort when cube-integral validation exceeds --validation_rtol.",
    )
    ap.add_argument(
        "--validation_rtol",
        type=float,
        default=0.05,
        help="Relative error threshold used with --fail_on_validation. Default is 0.05.",
    )
    return ap


def main(argv: Optional[List[str]] = None) -> None:
    args = build_parser().parse_args(argv)

    total_timer = Timer()
    eprint(f"[Start] {args.log.name}  state={args.state}")
    eprint(f"[Implementation] {IMPLEMENTATION_NOTE}")

    if not args.log.exists():
        raise FileNotFoundError(f"Log file not found: {args.log}")
    if args.chk is not None and not args.chk.exists():
        raise FileNotFoundError(f"Checkpoint file not found: {args.chk}")
    if not args.cpl_only:
        if args.state is None:
            raise ValueError("--state is required unless --cpl_only is specified.")
        if args.mocubes_dir is None and args.chk is None:
            raise ValueError("Provide either --chk or --mocubes_dir unless --cpl_only is specified.")

    dtype = np.float32 if args.dtype == "float32" else np.float64

    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    mode_tag = args.density_mode.replace("-", "_")
    if args.density_mode == "contribution-map" and not args.skip_validation_report:
        eprint("[Mode warning] contribution-map mode uses signed square 2*c*abs(c) weights to display configuration contributions. Moment-integral validation against Gaussian is not expected to pass in this mode.")
    if args.density_mode == "validation":
        eprint(
            "[Mode] validation: EDTM weight = "
            f"{args.validation_edtm_factor:g}*c; MDTM scale = {args.validation_mdtm_factor:g}*c; "
            f"phase_align={args.phase_align}"
        )
    else:
        eprint("[Mode] contribution-map: EDTM weight = 2*c*abs(c); MDTM scale = c*abs(c).")

    mode_settings_path = outdir / "density_mode_settings.txt"
    mode_settings_path.write_text(
        "density_mode=" + args.density_mode + "\n"
        + "phase_align=" + args.phase_align + "\n"
        + f"validation_edtm_factor={args.validation_edtm_factor}\n"
        + f"validation_mdtm_factor={args.validation_mdtm_factor}\n"
        + "accumulation=" + args.accumulation + "\n"
        + "implementation_origin=independent_python_from_published_equations\n"
        + "matlab_dependency=none\n"
        + "contribution_map_rationale=signed_square_weights_follow_configuration_contributions_not_linear_transition_moments\n",
        encoding="utf-8",
    )

    t0 = time.perf_counter()
    log_text = args.log.read_text(errors="ignore")

    cpl_rows: Optional[List[CPLRow]] = None
    cpl_csv_path = args.cpl_csv_path or (outdir / f"{args.log.stem}-CPL.csv")
    if not args.no_cpl_csv:
        cpl_rows = build_cpl_rows(log_text)
        write_cpl_csv(cpl_csv_path, args.log, cpl_rows)
        eprint(
            f"[1/8] Wrote CPL CSV in {Timer.fmt(time.perf_counter() - t0)}: "
            f"states={len(cpl_rows)} path={cpl_csv_path}"
        )
        t0 = time.perf_counter()

    if args.cpl_only:
        eprint("[Done] CPL CSV export only")
        eprint(f"  Total elapsed: {Timer.fmt(total_timer.elapsed())}")
        if not args.no_cpl_csv:
            eprint(f"  CPL CSV output: {cpl_csv_path}")
        return

    excitations = parse_excited_state_block(log_text, args.state)
    coeff_diag = diagnose_coefficients(log_text, args.state, excitations)
    eprint(
        "[CI] "
        f"state={coeff_diag.state} n={coeff_diag.n_configurations} "
        f"2*sum(c^2)={coeff_diag.two_sum_c_squared:.5f} "
        f"IOp(9/40)={coeff_diag.iop_940_level}"
    )
    if coeff_diag.warning:
        eprint(f"[CI warning] {coeff_diag.warning}")
    if args.strict_ci and (not coeff_diag.has_recommended_iop or coeff_diag.two_sum_c_squared < 0.98):
        raise ValueError(
            "The printed CI coefficients did not pass the strict diagnostic. "
            "Repeat the Gaussian calculation with IOp(9/40=4), then rerun CPLkit."
        )

    mux, muy, muz = parse_edtm_vector(log_text, args.state)
    mx, my, mz = parse_mdtm_vector(log_text, args.state)
    mu_norm = math.sqrt(mux * mux + muy * muy + muz * muz)
    m_norm = math.sqrt(mx * mx + my * my + mz * mz)
    if mu_norm == 0.0:
        raise ValueError(f"EDTM vector norm is zero for state {args.state}.")
    if m_norm == 0.0:
        raise ValueError(f"MDTM vector norm is zero for state {args.state}.")
    eprint(
        f"[2/8] Parsed target state in {Timer.fmt(time.perf_counter() - t0)}: "
        f"excitations={len(excitations)} |mu|={mu_norm:.6f} Au |m|={m_norm:.6f} Au"
    )

    outprefix = args.outprefix or f"S{args.state}"
    logstem = args.log.stem
    baseprefix = f"{logstem}_{outprefix}"
    density_baseprefix = f"{baseprefix}_{mode_tag}"
    needed_mos = sorted({ex.occ for ex in excitations} | {ex.virt for ex in excitations})
    eprint(f"[Info] Needed MOs: {needed_mos}")
    mo_cube_paths = ensure_mo_cubes(needed_mos, args, baseprefix, outdir)

    t0 = time.perf_counter()
    ref_cube = read_cube(mo_cube_paths[needed_mos[0]], dtype=dtype, cube_units=args.cube_units)
    eprint(
        f"[4/8] Read reference cube in {Timer.fmt(time.perf_counter() - t0)} "
        f"grid=({ref_cube.nx},{ref_cube.ny},{ref_cube.nz}) nvals={ref_cube.nvals} val_ids={ref_cube.val_ids}"
    )

    coords_mode = args.coords
    if coords_mode == "auto":
        try:
            _ = coordinate_1d_aligned_bohr(ref_cube)
            coords_mode = "aligned"
        except Exception:
            coords_mode = "general"
    eprint(f"[Info] coords_mode={coords_mode}")
    eprint(f"[Info] cube_units={args.cube_units}  unit_scale_to_bohr={ref_cube.unit_scale_to_bohr:g}")

    phi_cache: "OrderedDict[int, np.ndarray]" = OrderedDict()

    def get_phi(mo: int) -> np.ndarray:
        if mo in phi_cache:
            phi_cache.move_to_end(mo)
            return phi_cache[mo]

        c = read_cube(mo_cube_paths[mo], dtype=dtype, cube_units=args.cube_units)
        if (c.nx, c.ny, c.nz) != (ref_cube.nx, ref_cube.ny, ref_cube.nz):
            raise ValueError(f"Grid mismatch for mo{mo}")
        if not np.isclose(c.unit_scale_to_bohr, ref_cube.unit_scale_to_bohr):
            raise ValueError(f"Unit conversion mismatch for mo{mo}. Regenerate cubes consistently.")
        if c.nvals != ref_cube.nvals:
            raise ValueError(f"nvals mismatch for mo{mo}. Regenerate cubes consistently.")

        phi_cache[mo] = c.data.astype(dtype, copy=False)
        phi_cache.move_to_end(mo)
        while args.max_phi_cache >= 0 and len(phi_cache) > args.max_phi_cache:
            phi_cache.popitem(last=False)
        return phi_cache[mo]

    rho_edtm_x = np.zeros((ref_cube.nx, ref_cube.ny, ref_cube.nz), dtype=dtype)
    rho_edtm_y = np.zeros_like(rho_edtm_x)
    rho_edtm_z = np.zeros_like(rho_edtm_x)

    rho_mdtm_x = np.zeros_like(rho_edtm_x)
    rho_mdtm_y = np.zeros_like(rho_edtm_x)
    rho_mdtm_z = np.zeros_like(rho_edtm_x)

    eprint("[5/8] Accumulating excitations for EDTM and MDTM ...")
    eprint(f"[Info] requested_accumulation={args.accumulation}  vectorized_x_block={args.vectorized_x_block}  stack_mode={args.vectorized_stack_mode}")
    ex_t0 = time.perf_counter()
    n_ex = len(excitations)

    if args.accumulation == "vectorized":
        if coords_mode != "aligned":
            msg = "[Vectorized] Non-aligned grid detected; vectorized accumulation cannot be used."
            if args.require_vectorized:
                raise ValueError(msg + " Regenerate cubes on an axis-aligned grid or run with --accumulation stream.")
            eprint(msg + " Falling back to streaming accumulation.")
            args.accumulation = "stream"
        else:
            eprint("[Vectorized] Enabled. Progress should be reported as 'stack' and 'blocks', not as 'excit'.")
            scratch_dir = args.vectorized_scratch_dir or (outdir / "_scratch")
            stack_label = f"{baseprefix}_stack_{args.dtype}_{len(needed_mos)}mos.dat"
            accumulate_vectorized_blocks(
                ref_cube,
                coords_mode,
                excitations,
                get_phi,
                dtype,
                (rho_edtm_x, rho_edtm_y, rho_edtm_z),
                (rho_mdtm_x, rho_mdtm_y, rho_mdtm_z),
                x_block=args.vectorized_x_block,
                stack_mode=args.vectorized_stack_mode,
                scratch_dir=scratch_dir,
                stack_label=stack_label,
                progress=eprint,
                density_mode=args.density_mode,
                validation_edtm_factor=args.validation_edtm_factor,
                validation_mdtm_factor=args.validation_mdtm_factor,
            )

    if args.accumulation == "stream":
        for i, ex in enumerate(excitations, 1):
            if i == 1 or i % 5 == 0 or i == n_ex:
                eprint(progress_line("  excit", i - 1, n_ex, ex_t0))

            w_edtm, s = excitation_weights(
                ex.coeff,
                args.density_mode,
                validation_edtm_factor=args.validation_edtm_factor,
                validation_mdtm_factor=args.validation_mdtm_factor,
            )
            phi_occ = get_phi(ex.occ)
            phi_virt = get_phi(ex.virt)

            accumulate_edtm_components(
                ref_cube,
                coords_mode,
                phi_occ,
                phi_virt,
                w_edtm,
                (rho_edtm_x, rho_edtm_y, rho_edtm_z),
            )

            grad_occ = orbital_gradients(phi_occ, ref_cube)
            grad_virt = orbital_gradients(phi_virt, ref_cube)
            accumulate_mdtm_components(
                ref_cube,
                coords_mode,
                phi_occ,
                phi_virt,
                grad_occ,
                grad_virt,
                s,
                (rho_mdtm_x, rho_mdtm_y, rho_mdtm_z),
            )

        eprint(progress_line("  excit", n_ex, n_ex, ex_t0))
    else:
        eprint(f"[Vectorized] Finished block accumulation in {Timer.fmt(time.perf_counter() - ex_t0)}")

    if args.phase_align == "edtm":
        vol = cube_voxel_volume_bohr3(ref_cube)
        edtm_integral = (
            float(np.sum(rho_edtm_x, dtype=np.float64) * vol),
            float(np.sum(rho_edtm_y, dtype=np.float64) * vol),
            float(np.sum(rho_edtm_z, dtype=np.float64) * vol),
        )
        dot = edtm_integral[0] * mux + edtm_integral[1] * muy + edtm_integral[2] * muz
        eprint(f"[Phase] EDTM integral before phase alignment = {edtm_integral}; dot_with_Gaussian={dot:.8g}")
        if dot < 0.0:
            eprint("[Phase] Multiplying all reconstructed EDTM/MDTM densities by -1 to align EDTM with Gaussian.")
            rho_edtm_x *= -1
            rho_edtm_y *= -1
            rho_edtm_z *= -1
            rho_mdtm_x *= -1
            rho_mdtm_y *= -1
            rho_mdtm_z *= -1

    validation_rows = []
    if not args.skip_validation_report:
        validation_rows.extend(
            compare_component_integrals(ref_cube, "EDTM", (rho_edtm_x, rho_edtm_y, rho_edtm_z), (mux, muy, muz))
        )
        validation_rows.extend(
            compare_component_integrals(ref_cube, "MDTM", (rho_mdtm_x, rho_mdtm_y, rho_mdtm_z), (mx, my, mz))
        )
        validation_csv_path = args.validation_csv_path or (outdir / f"{density_baseprefix}_moment_validation.csv")
        validation_csv_path.parent.mkdir(parents=True, exist_ok=True)
        write_validation_csv(validation_csv_path, validation_rows)
        eprint(f"[Validation] Wrote moment-density validation: {validation_csv_path}")
        worst = max((row.relative_error for row in validation_rows if math.isfinite(row.relative_error)), default=0.0)
        if args.fail_on_validation and worst > args.validation_rtol:
            raise ValueError(
                f"Moment-density validation failed: max relative error {worst:.4g} exceeds {args.validation_rtol:.4g}."
            )

    eprint("[6/8] Projecting onto transition vectors ...")
    muxn, muyn, muzn = mux / mu_norm, muy / mu_norm, muz / mu_norm
    mxn, myn, mzn = mx / m_norm, my / m_norm, mz / m_norm

    rho_edtm_total = muxn * rho_edtm_x + muyn * rho_edtm_y + muzn * rho_edtm_z
    rho_mdtm_total = mxn * rho_mdtm_x + myn * rho_mdtm_y + mzn * rho_mdtm_z

    eprint("[7/8] Writing cube files ...")
    edtm_comment1 = f"EDTM density cube from {args.log.name} (state {args.state}; mode={args.density_mode})"
    edtm_comment2 = f"Projected EDTM density; dtype={args.dtype}; mode={args.density_mode}"
    mdtm_comment1 = f"MDTM density cube from {args.log.name} (state {args.state}; mode={args.density_mode})"
    mdtm_comment2 = f"Projected MDTM density; dtype={args.dtype}; mode={args.density_mode}"

    edtm_total_path = outdir / f"{density_baseprefix}_EDTM_total.cube"
    mdtm_total_path = outdir / f"{density_baseprefix}_MDTM_total.cube"

    write_cube(edtm_total_path, ref_cube, rho_edtm_total, edtm_comment1, edtm_comment2)
    write_cube(mdtm_total_path, ref_cube, rho_mdtm_total, mdtm_comment1, mdtm_comment2)

    if args.keep_components:
        write_cube(outdir / f"{density_baseprefix}_EDTM_x.cube", ref_cube, rho_edtm_x, edtm_comment1, edtm_comment2)
        write_cube(outdir / f"{density_baseprefix}_EDTM_y.cube", ref_cube, rho_edtm_y, edtm_comment1, edtm_comment2)
        write_cube(outdir / f"{density_baseprefix}_EDTM_z.cube", ref_cube, rho_edtm_z, edtm_comment1, edtm_comment2)

        write_cube(outdir / f"{density_baseprefix}_MDTM_x.cube", ref_cube, rho_mdtm_x, mdtm_comment1, mdtm_comment2)
        write_cube(outdir / f"{density_baseprefix}_MDTM_y.cube", ref_cube, rho_mdtm_y, mdtm_comment1, mdtm_comment2)
        write_cube(outdir / f"{density_baseprefix}_MDTM_z.cube", ref_cube, rho_mdtm_z, mdtm_comment1, mdtm_comment2)

    eprint("[8/8] Done")
    eprint(f"  Total elapsed: {Timer.fmt(total_timer.elapsed())}")
    eprint(f"  EDTM output: {edtm_total_path}")
    eprint(f"  MDTM output: {mdtm_total_path}")
    if not args.no_cpl_csv:
        eprint(f"  CPL CSV output: {cpl_csv_path}")


if __name__ == "__main__":
    main()
