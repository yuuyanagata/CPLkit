"""Gaussian cubegen discovery and execution helpers."""

from __future__ import annotations

import argparse
import hashlib
import os
import re
import shlex
import shutil
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Sequence, Tuple, Optional

import numpy as np

from .utils import eprint, progress_line


def candidate_cubegen_paths() -> List[str]:
    candidates: List[str] = []

    for name in ("cubegen", "cubegen.exe"):
        found = shutil.which(name)
        if found:
            candidates.append(found)

    env_vars = ["GAUSS_EXEDIR", "g16root", "g09root"]
    for var in env_vars:
        val = os.environ.get(var)
        if not val:
            continue
        parts = [p for p in re.split(r"[;:]", val) if p.strip()]
        for base in parts:
            p = Path(base)
            if p.is_file() and p.name.lower().startswith("cubegen"):
                candidates.append(str(p))
                continue
            for sub in (p, p / "g16", p / "g09"):
                for exe in ("cubegen", "cubegen.exe"):
                    cand = sub / exe
                    if cand.exists():
                        candidates.append(str(cand))

    common_windows = [
        Path(r"C:\Gaussian\g16\cubegen.exe"),
        Path(r"C:\G16W\g16\cubegen.exe"),
        Path(r"C:\Program Files\Gaussian\g16\cubegen.exe"),
        Path(r"C:\Program Files\g16\cubegen.exe"),
        Path(r"C:\Program Files\Gaussian 16W\g16\cubegen.exe"),
    ]
    for p in common_windows:
        if p.exists():
            candidates.append(str(p))

    unique: List[str] = []
    seen = set()
    for c in candidates:
        key = str(Path(c)).lower()
        if key not in seen:
            unique.append(c)
            seen.add(key)
    return unique


def resolve_cubegen(user_value: str) -> str:
    p = Path(user_value)
    if p.exists():
        return str(p)

    found = shutil.which(user_value)
    if found:
        return found

    if user_value in {"cubegen", "cubegen.exe"}:
        candidates = candidate_cubegen_paths()
        if candidates:
            return candidates[0]

    msg = [f"Could not find cubegen executable: {user_value}"]
    cand = candidate_cubegen_paths()
    if cand:
        msg.append("Candidate paths found:")
        msg.extend(f"  - {x}" for x in cand)
        msg.append("Use --cubegen with one of the paths above.")
    else:
        msg.append("No cubegen candidate paths were found.")
        msg.append(r'Please specify the full path with --cubegen "C:\Gaussian\g16\cubegen.exe"')
    raise FileNotFoundError("\n".join(msg))


def build_cubegen_command(
    cubegen_exe: str,
    npts: int,
    mo: int,
    chk_path: Path,
    cube_path: Path,
    grid_args: Sequence[str],
    prefix_args: Sequence[str] = (),
) -> List[str]:
    """Build a cubegen-compatible command.

    Gaussian cubegen accepts the form
        cubegen nprocs property fchk cube npts format

    Gaussian-compatible cube generator accepts the same positional arguments, but its own options such as
    --backend, --threads, and --grid must be placed before the nprocs argument.
    The optional prefix_args parameter is therefore inserted immediately after the
    executable path.
    """
    cmd = [cubegen_exe]
    cmd.extend(prefix_args)
    cmd.extend([str(npts), f"MO={mo}", str(chk_path), str(cube_path)])
    cmd.extend(grid_args)
    return cmd




def build_multi_mo_cubegen_command(
    cubegen_exe: str,
    npts: int,
    mos: Sequence[int],
    chk_path: Path,
    cube_path: Path,
    grid_args: Sequence[str],
    prefix_args: Sequence[str] = (),
) -> List[str]:
    """Build a Gaussian-compatible cube generator-compatible multi-MO cube-generation command."""
    if not mos:
        raise ValueError("mos must not be empty")
    mo_spec = ",".join(str(int(mo)) for mo in mos)
    cmd = [cubegen_exe]
    cmd.extend(prefix_args)
    cmd.extend([str(npts), f"MO={mo_spec}", str(chk_path), str(cube_path)])
    cmd.extend(grid_args)
    return cmd


def _replace_backend_for_fallback(cmd: Sequence[str], backend: str = "cpu") -> Optional[List[str]]:
    """Return a CPU-fallback command for Gaussian-compatible cube generator-style options.

    Gaussian-compatible cube generatorCpp can fail in CUDA mode for a particular MO or CUDA runtime
    condition.  The positional cubegen-compatible arguments remain valid for
    the CPU backend, so CPLkit can preserve the same grid and output path while
    changing only the backend selector.
    """
    cmd2 = list(cmd)
    for i, token in enumerate(cmd2):
        low = token.lower()
        if low == "--backend" and i + 1 < len(cmd2):
            current = cmd2[i + 1].lower()
            if current in {"cuda", "auto"}:
                cmd2[i + 1] = backend
                return cmd2
        if low.startswith("--backend="):
            current = token.split("=", 1)[1].lower()
            if current in {"cuda", "auto"}:
                cmd2[i] = "--backend=" + backend
                return cmd2
    return None


def _find_output_cube_from_command(cmd: Sequence[str]) -> Optional[Path]:
    """Best-effort extraction of the cube output path from a cubegen-compatible command."""
    for i, token in enumerate(cmd):
        if token.upper().startswith(("MO=", "AMO=", "BMO=", "HOMO", "LUMO", "DENSITY", "FDENSITY")):
            # executable, optional prefix args, nprocs, property, fchk, cube, grid...
            if i + 2 < len(cmd):
                return Path(cmd[i + 2])
    return None


def _format_cubegen_failure(cmd: Sequence[str], res: subprocess.CompletedProcess[str], *, prefix: str = "cubegen failed") -> str:
    msg = [
        f"{prefix} with exit status {res.returncode}.",
        "Command:",
        "  " + subprocess.list2cmdline(list(cmd)),
    ]
    if res.stdout and res.stdout.strip():
        msg.append("stdout:\n" + res.stdout.rstrip())
    if res.stderr and res.stderr.strip():
        msg.append("stderr:\n" + res.stderr.rstrip())
    if os.name == "nt" and res.returncode in {127, 3228369022}:
        msg.append(
            "On Windows, this exit status can indicate a missing runtime DLL, a CUDA runtime failure, "
            "or an unhandled native exception inside the cubegen-compatible executable."
        )
    return "\n".join(msg)


def run_cubegen(cmd: Sequence[str]) -> None:
    res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if res.returncode == 0:
        return

    fallback_cmd = _replace_backend_for_fallback(cmd, backend="cpu")
    if fallback_cmd is not None:
        cube_path = _find_output_cube_from_command(fallback_cmd)
        if cube_path is not None and cube_path.exists():
            try:
                cube_path.unlink()
            except OSError:
                pass
        eprint(
            "[WARNING] cubegen-compatible command failed in CUDA/auto mode. "
            "Retrying the same MO cube with CPU backend.\n"
            + _format_cubegen_failure(cmd, res, prefix="Initial cubegen attempt failed")
        )
        res2 = subprocess.run(fallback_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if res2.returncode == 0:
            eprint("[Info] CPU fallback for this MO cube completed successfully.")
            return
        raise RuntimeError(
            _format_cubegen_failure(cmd, res, prefix="Initial cubegen attempt failed")
            + "\n\n"
            + _format_cubegen_failure(fallback_cmd, res2, prefix="CPU fallback cubegen attempt also failed")
        )

    raise RuntimeError(_format_cubegen_failure(cmd, res))





def _parse_multi_cube_structure(path: Path) -> tuple[list[str], int, list[int], int]:
    """Parse a Gaussian-compatible cube generator/Gaussian multi-orbital cube header.

    The common format uses a negative atom count and an additional line after
    the atom block such as "2 127 128".  Some variants also include the number
    of values as the fifth field of the atom header.  Both conventions are
    accepted.
    """
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    if len(lines) < 7:
        raise ValueError(f"{path} does not look like a multi-orbital cube file.")
    atom_tokens = lines[2].split()
    if len(atom_tokens) < 4:
        raise ValueError(f"Bad cube atom header line in {path}: {lines[2]!r}")
    natoms_raw = int(atom_tokens[0])
    natoms = abs(natoms_raw)
    cursor = 6 + natoms
    if len(lines) <= cursor:
        raise ValueError(f"{path}: missing multi-orbital dataset line or cube data.")

    nvals_from_header = int(atom_tokens[4]) if len(atom_tokens) >= 5 else None
    dataset_tokens = lines[cursor].split()
    ids: list[int] = []

    if nvals_from_header is not None:
        nvals = nvals_from_header
        maybe = [int(x) for x in dataset_tokens] if dataset_tokens else []
        if len(maybe) == nvals + 1 and maybe[0] == nvals:
            ids = maybe[1:]
        elif len(maybe) >= nvals:
            ids = maybe[:nvals]
        else:
            ids = maybe
        cursor += 1
    else:
        maybe = [int(x) for x in dataset_tokens]
        if not maybe:
            raise ValueError(f"{path}: empty multi-orbital dataset line.")
        nvals = maybe[0]
        ids = maybe[1:]
        cursor += 1
        while len(ids) < nvals and cursor < len(lines):
            more = lines[cursor].split()
            if not more:
                cursor += 1
                continue
            try:
                ids.extend(int(x) for x in more)
                cursor += 1
            except ValueError:
                break

    if nvals < 1:
        raise ValueError(f"{path}: invalid number of cube values: {nvals}")
    if len(ids) < nvals:
        ids.extend(range(len(ids) + 1, nvals + 1))
    ids = ids[:nvals]

    # The split files must be ordinary scalar cube files, so the atom count is
    # made positive and the extra multi-MO dataset line is omitted.
    one_atom_line = f"{natoms:5d} {float(atom_tokens[1]): .6f} {float(atom_tokens[2]): .6f} {float(atom_tokens[3]): .6f}"
    one_header = [lines[0], lines[1], one_atom_line, lines[3], lines[4], lines[5]]
    one_header.extend(lines[6:6 + natoms])
    return one_header, nvals, ids, cursor


def split_multi_mo_cube(multi_cube_path: Path, mo_to_path: Dict[int, Path], dtype_name: str = "float64") -> Dict[int, Path]:
    """Split a multi-orbital cube into ordinary mo<N>.cube files.

    Gaussian-compatible cube generatorCpp writes multi-orbital values point-major, for example
    MO127(point1) MO128(point1) MO127(point2) MO128(point2) ... .
    CPLkit keeps scalar mo<N>.cube files as its on-disk cache format, so this
    splitter preserves compatibility with the rest of the analysis pipeline.
    """
    one_header, nvals, ids, data_cursor = _parse_multi_cube_structure(multi_cube_path)
    lines = multi_cube_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    nx = abs(int(one_header[3].split()[0]))
    ny = abs(int(one_header[4].split()[0]))
    nz = abs(int(one_header[5].split()[0]))
    npoints = nx * ny * nz
    expected = npoints * nvals
    dtype = np.float32 if str(dtype_name).lower() == "float32" else np.float64
    data_text = "\n".join(lines[data_cursor:]).strip()
    arr = np.fromstring(data_text, sep=" ", dtype=dtype)
    if arr.size < expected:
        raise ValueError(f"{multi_cube_path}: cube data short. expected {expected}, got {arr.size}.")
    if arr.size > expected:
        arr = arr[:expected]
    values = arr.reshape((npoints, nvals), order="C")
    id_to_col = {int(mo_id): idx for idx, mo_id in enumerate(ids)}
    written: Dict[int, Path] = {}
    for mo, out_path in mo_to_path.items():
        if int(mo) not in id_to_col:
            raise ValueError(f"{multi_cube_path}: requested MO {mo} not found in multi-cube dataset IDs {ids}")
        col = values[:, id_to_col[int(mo)]]
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with out_path.open("w", encoding="utf-8") as f:
            f.write(one_header[0] + "\n")
            f.write(f"Split from {multi_cube_path.name}; MO={mo}\n")
            for line in one_header[2:]:
                f.write(line.rstrip("\n") + "\n")
            for i in range(0, npoints, 6):
                chunk = col[i:i + 6]
                f.write(" ".join(f"{float(v): .5E}" for v in chunk) + "\n")
        written[int(mo)] = out_path
    return written


def _safe_token(text: str, max_len: int = 96) -> str:
    token = re.sub(r"[^A-Za-z0-9]+", "_", text).strip("_")
    return (token[:max_len].strip("_") or "default").lower()


def _settings_token(cubegen_exe: str, prefix_args: Sequence[str], npts: int, grid_args: Sequence[str]) -> str:
    """Return a short but collision-resistant token for an MO-cube directory.

    Earlier validation packages embedded the full cubegen command in the
    directory name.  That protected against mixing cubes from different grids,
    but it could exceed the classic Windows MAX_PATH limit when the package was
    unpacked below a long working directory.  This compact token preserves the
    relevant human-readable settings and appends a SHA1 digest of the full
    cubegen command settings.
    """
    engine_raw = Path(cubegen_exe).stem.lower()
    if "multi-mo-cubegen" in engine_raw or "multimocubegen" in engine_raw:
        engine = "fcpp"
    elif "multimocubegen" in engine_raw:
        engine = "fcg"
    elif "cubegen" in engine_raw:
        engine = "cubegen"
    else:
        engine = _safe_token(engine_raw, 12)

    raw = "__".join([
        str(Path(cubegen_exe)),
        "prefix", " ".join(prefix_args) or "none",
        "nprocs", str(npts),
        "grid", " ".join(grid_args) or "none",
    ])
    digest = hashlib.sha1(raw.encode("utf-8", errors="ignore")).hexdigest()[:10]

    grid_part = _safe_token("_".join(grid_args) or "default", 18)
    backend = "auto"
    threads = None
    for i, token in enumerate(prefix_args):
        token_l = token.lower()
        if token_l == "--backend" and i + 1 < len(prefix_args):
            backend = _safe_token(prefix_args[i + 1], 8)
        elif token_l.startswith("--backend="):
            backend = _safe_token(token.split("=", 1)[1], 8)
        elif token_l == "--threads" and i + 1 < len(prefix_args):
            threads = _safe_token(prefix_args[i + 1], 6)
        elif token_l.startswith("--threads="):
            threads = _safe_token(token.split("=", 1)[1], 6)

    parts = [engine, f"b{backend}"]
    if threads:
        parts.append(f"t{threads}")
    parts.extend([f"n{npts}", f"g{grid_part}", digest])
    return "_".join(parts)


def _cube_header_signature(path: Path, ndigits: int = 10) -> Tuple:
    """Return a compact cube-header signature for grid compatibility checks.

    The signature intentionally ignores the two comment lines because Gaussian cubegen
    and Gaussian-compatible cube generator may write different comments for numerically compatible grids.
    """
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        lines = []
        for _ in range(256):
            line = f.readline()
            if line == "":
                break
            lines.append(line.rstrip("\n"))
            if len(lines) >= 6:
                try:
                    natoms_raw = int(lines[2].split()[0])
                    natoms = abs(natoms_raw)
                    needed = 6 + natoms
                    if len(lines) >= needed:
                        # Optional dataset-id line is not needed for grid compatibility.
                        break
                except Exception:
                    pass
    if len(lines) < 6:
        raise ValueError(f"{path} does not look like a cube file; fewer than six header lines were found.")

    def parse_header(i: int):
        toks = lines[i].split()
        return (int(toks[0]), tuple(round(float(x), ndigits) for x in toks[1:4]))

    atom = parse_header(2)
    x = parse_header(3)
    y = parse_header(4)
    z = parse_header(5)
    natoms = abs(atom[0])
    atom_lines = tuple(lines[6:6 + natoms])
    nvals = None
    toks = lines[2].split()
    if len(toks) >= 5:
        nvals = int(toks[4])
    return (atom, x, y, z, nvals, atom_lines)


def _group_by_cube_header(paths: Dict[int, Path]) -> Dict[Tuple, List[int]]:
    groups: Dict[Tuple, List[int]] = {}
    for mo, path in paths.items():
        sig = _cube_header_signature(path)
        groups.setdefault(sig, []).append(mo)
    return groups


def _format_mo_list(mos: List[int], limit: int = 24) -> str:
    mos = sorted(mos)
    shown = ", ".join(str(x) for x in mos[:limit])
    if len(mos) > limit:
        shown += f", ... ({len(mos)} total)"
    return shown


def _validate_or_repair_cube_headers(
    mo_cube_paths: Dict[int, Path],
    *,
    generated_mos: List[int],
    args: argparse.Namespace,
    pending_builder,
    runner,
) -> Dict[int, Path]:
    """Ensure that all MO cubes share an identical grid header.

    If CPLkit controls the cube directory and a checkpoint is available, a mixed
    existing directory is repaired by regenerating the inconsistent cubes.  If a
    user supplied --mocubes_dir, the function raises an explicit error because
    regenerating user-managed cubes would be surprising.
    """
    if len(mo_cube_paths) <= 1:
        return mo_cube_paths

    groups = _group_by_cube_header(mo_cube_paths)
    if len(groups) == 1:
        return mo_cube_paths

    eprint("[WARNING] Existing MO cubes use inconsistent cube grids.")
    for idx, (_, mos) in enumerate(sorted(groups.items(), key=lambda kv: len(kv[1]), reverse=True), 1):
        eprint(f"  grid group {idx}: MOs {_format_mo_list(mos)}")

    if args.mocubes_dir is not None:
        raise ValueError(
            "Grid mismatch was detected in --mocubes_dir. Regenerate the MO cubes in a clean folder, "
            "or provide a folder whose mo<N>.cube files were generated with one identical grid setting."
        )

    if args.chk is None:
        raise ValueError("Grid mismatch was detected, but no checkpoint file was supplied for regeneration.")

    desired_sig: Optional[Tuple] = None
    for mo in generated_mos:
        if mo in mo_cube_paths:
            desired_sig = _cube_header_signature(mo_cube_paths[mo])
            break

    if desired_sig is None:
        # All cubes were pre-existing. Because the current cubegen settings are known,
        # the least ambiguous repair is to regenerate all required cubes.
        eprint("[WARNING] All cubes were pre-existing and mixed. Regenerating all required MO cubes with the current settings.")
        repair_mos = sorted(mo_cube_paths)
    else:
        repair_mos = sorted(mo for mo, path in mo_cube_paths.items() if _cube_header_signature(path) != desired_sig)
        eprint(f"[WARNING] Regenerating inconsistent existing MO cubes: {_format_mo_list(repair_mos)}")

    for mo in repair_mos:
        try:
            mo_cube_paths[mo].unlink()
        except FileNotFoundError:
            pass

    items = [pending_builder(mo, mo_cube_paths[mo]) for mo in repair_mos]
    repaired = runner(items, label="  cubegen-repair")
    mo_cube_paths.update(repaired)

    groups = _group_by_cube_header(mo_cube_paths)
    if len(groups) != 1:
        raise ValueError(
            "Grid mismatch remains after regeneration. Delete the relevant *_mocubes folder and rerun, "
            "or set OVERWRITE_MO_CUBES=1 for one clean run."
        )
    return mo_cube_paths

def ensure_mo_cubes(
    needed_mos: List[int],
    args: argparse.Namespace,
    baseprefix: str,
    outdir: Path,
) -> Dict[int, Path]:
    mo_cube_paths: Dict[int, Path] = {}

    if args.mocubes_dir is not None:
        if not args.mocubes_dir.exists():
            raise FileNotFoundError(f"MO cube directory not found: {args.mocubes_dir}")
        for mo in needed_mos:
            p = args.mocubes_dir / f"mo{mo}.cube"
            if not p.exists():
                raise FileNotFoundError(f"Missing MO cube: {p}")
            mo_cube_paths[mo] = p
        eprint(f"[3/8] Using existing MO cubes in {args.mocubes_dir}")
        groups = _group_by_cube_header(mo_cube_paths)
        if len(groups) != 1:
            eprint("[WARNING] The supplied --mocubes_dir contains mixed cube grids.")
            for idx, (_, mos) in enumerate(sorted(groups.items(), key=lambda kv: len(kv[1]), reverse=True), 1):
                eprint(f"  grid group {idx}: MOs {_format_mo_list(mos)}")
            raise ValueError(
                "Grid mismatch was detected in --mocubes_dir. Use a clean MO cube folder generated with one grid setting."
            )
        return mo_cube_paths

    if args.chk is None:
        raise ValueError("Provide either --chk or --mocubes_dir.")

    cubegen_exe = resolve_cubegen(args.cubegen)
    grid_args = shlex.split(args.cubegen_grid, posix=False)
    prefix_args = shlex.split(getattr(args, "cubegen_prefix", "") or "", posix=False)
    settings_token = _settings_token(cubegen_exe, prefix_args, args.cubegen_npts, grid_args)
    base_token = _safe_token(baseprefix, 40)
    mocdir = outdir / f"{base_token}_mc_{settings_token}"
    try:
        mocdir.mkdir(parents=True, exist_ok=True)
    except FileNotFoundError as exc:
        # On Windows this often indicates that the fully qualified path exceeded
        # the legacy MAX_PATH limit.  Fall back to an even shorter directory name.
        short_base = hashlib.sha1(baseprefix.encode("utf-8", errors="ignore")).hexdigest()[:8]
        mocdir = outdir / f"mc_{short_base}_{settings_token.split('_')[-1]}"
        mocdir.mkdir(parents=True, exist_ok=True)

    jobs = max(1, int(getattr(args, "cubegen_jobs", 1) or 1))
    multi_chunk = max(1, int(getattr(args, "cubegen_multi_mo_chunk", 1) or 1))

    mocdir.mkdir(parents=True, exist_ok=True)
    meta_path = mocdir / "mocubes_settings.txt"
    if not meta_path.exists() or args.overwrite_mo_cubes:
        meta_path.write_text(
            "baseprefix=" + baseprefix + "\n"
            + "cubegen=" + cubegen_exe + "\n"
            + "cubegen_prefix=" + " ".join(prefix_args) + "\n"
            + "cubegen_nprocs=" + str(args.cubegen_npts) + "\n"
            + "cubegen_grid=" + args.cubegen_grid + "\n"
            + "cubegen_multi_mo_chunk=" + str(multi_chunk) + "\n"
            + "multi_mo_output_mode=separate_per_orbital\n"
            + "settings_token=" + settings_token + "\n",
            encoding="utf-8",
        )

    eprint("[3/8] Generating MO cubes with cubegen ...")
    eprint(f"[Info] cubegen = {cubegen_exe}")
    eprint(f"[Info] cubegen_prefix = {' '.join(prefix_args) if prefix_args else '(none)'}")
    eprint(f"[Info] cubegen_nprocs = {args.cubegen_npts}")
    eprint(f"[Info] cubegen_grid = {args.cubegen_grid}")
    eprint(f"[Info] cubegen_jobs = {jobs}")
    eprint(f"[Info] cubegen_multi_mo_chunk = {multi_chunk}")
    eprint(f"[Info] mocubes_dir = {mocdir}")
    gen_t0 = time.perf_counter()

    pending: List[tuple[int, Path, List[str]]] = []
    for mo in needed_mos:
        cube_path = mocdir / f"mo{mo}.cube"
        if cube_path.exists() and not args.overwrite_mo_cubes:
            mo_cube_paths[mo] = cube_path
            continue
        cmd = build_cubegen_command(
            cubegen_exe=cubegen_exe,
            npts=args.cubegen_npts,
            mo=mo,
            chk_path=args.chk,
            cube_path=cube_path,
            grid_args=grid_args,
            prefix_args=prefix_args,
        )
        pending.append((mo, cube_path, cmd))

    if not pending:
        eprint("  cubegen: all required MO cubes already exist; checking header consistency before reuse")
    else:
        if multi_chunk > 1:
            nchunks = (len(pending) + multi_chunk - 1) // multi_chunk
            worker_count = min(jobs, nchunks)
            eprint(f"[Info] cubegen_pending_mos = {len(pending)}")
            eprint(f"[Info] cubegen_pending_multi_mo_chunks = {nchunks}")
            eprint(f"[Info] cubegen_max_workers = {worker_count}")
            eprint(
                f"[Info] Multi-MO cubegen is enabled: each multi-MO cubegen process will request up to "
                f"{multi_chunk} orbitals and Gaussian-compatible cube generatorCpp will write one cube file per orbital."
            )
        else:
            worker_count = min(jobs, len(pending))
            eprint(f"[Info] cubegen_pending = {len(pending)}")
            eprint(f"[Info] cubegen_max_workers = {worker_count}")
        if worker_count > 1:
            eprint(f"[Info] Parallel cubegen is enabled: up to {worker_count} cubegen subprocesses will be active concurrently.")
        else:
            eprint("[Info] Parallel cubegen is disabled because cubegen_jobs=1 or only one cube/chunk is pending.")

    def _build_item(mo: int, cube_path: Path) -> tuple[int, Path, List[str]]:
        cmd = build_cubegen_command(
            cubegen_exe=cubegen_exe,
            npts=args.cubegen_npts,
            mo=mo,
            chk_path=args.chk,
            cube_path=cube_path,
            grid_args=grid_args,
            prefix_args=prefix_args,
        )
        return (mo, cube_path, cmd)

    def _run_one(item: tuple[int, Path, List[str]]) -> tuple[int, Path]:
        mo, cube_path, cmd = item
        run_cubegen(cmd)
        if not cube_path.exists():
            raise FileNotFoundError(f"cubegen reported success but output cube was not created: {cube_path}")
        return mo, cube_path

    def _run_many(items: List[tuple[int, Path, List[str]]], label: str = "  cubegen") -> Dict[int, Path]:
        if not items:
            return {}
        local_t0 = time.perf_counter()
        local_workers = min(jobs, len(items))
        results: Dict[int, Path] = {}
        if jobs == 1 or len(items) == 1:
            for idx, item in enumerate(items, 1):
                eprint(progress_line(label, idx - 1, len(items), local_t0))
                mo, cube_path = _run_one(item)
                results[mo] = cube_path
        else:
            completed = 0
            eprint(progress_line(label, completed, len(items), local_t0))
            with ThreadPoolExecutor(max_workers=local_workers) as executor:
                future_to_mo = {executor.submit(_run_one, item): item[0] for item in items}
                for future in as_completed(future_to_mo):
                    mo, cube_path = future.result()
                    results[mo] = cube_path
                    completed += 1
                    eprint(progress_line(label, completed, len(items), local_t0))
        eprint(progress_line(label, len(items), len(items), local_t0))
        return results

    def _chunk_pending(items: List[tuple[int, Path, List[str]]]) -> List[tuple[List[int], Dict[int, Path], List[str], Path]]:
        chunks: List[tuple[List[int], Dict[int, Path], List[str], Path]] = []
        for chunk_index in range(0, len(items), multi_chunk):
            chunk = items[chunk_index:chunk_index + multi_chunk]
            mos = [mo for mo, _, _ in chunk]
            paths = {mo: cube_path for mo, cube_path, _ in chunk}
            digest = hashlib.sha1(",".join(str(x) for x in mos).encode("ascii")).hexdigest()[:10]
            tmp_path = mocdir / f"multi_{mos[0]}_{mos[-1]}_{len(mos)}_{digest}.cube"
            if tmp_path.exists():
                try:
                    tmp_path.unlink()
                except OSError:
                    pass
            cmd = build_multi_mo_cubegen_command(
                cubegen_exe=cubegen_exe,
                npts=args.cubegen_npts,
                mos=mos,
                chk_path=args.chk,
                cube_path=tmp_path,
                grid_args=grid_args,
                prefix_args=prefix_args,
            )
            chunks.append((mos, paths, cmd, tmp_path))
        return chunks

    def _expected_separate_multi_outputs(base_cube_path: Path, mos: Sequence[int]) -> Dict[int, Path]:
        # The revised Gaussian-compatible cube generatorCpp writes one scalar cube per requested orbital.
        # Example: frontier.cube with MO=127,128 produces frontier_MO127.cube and frontier_MO128.cube.
        return {int(mo): base_cube_path.with_name(f"{base_cube_path.stem}_MO{int(mo)}{base_cube_path.suffix}") for mo in mos}

    def _run_one_multi(item: tuple[List[int], Dict[int, Path], List[str], Path]) -> Dict[int, Path]:
        mos, paths, cmd, tmp_path = item
        expected = _expected_separate_multi_outputs(tmp_path, mos)

        # Remove stale outputs before launching Gaussian-compatible cube generator.  For multi-orbital
        # chunks, the revised MIT Gaussian-compatible cube generator writes one file per orbital, for
        # example:
        #     multi_127_134_8_xxxxx_MO127.cube
        # For a chunk containing only one MO, however, the executable may write
        # the requested base output file itself:
        #     multi_146_146_1_xxxxx.cube
        # rather than the suffixed form:
        #     multi_146_146_1_xxxxx_MO146.cube
        # CPLkit therefore accepts both forms for single-MO chunks.
        cleanup_paths = set(expected.values())
        cleanup_paths.add(tmp_path)
        for produced_path in cleanup_paths:
            if produced_path.exists():
                try:
                    produced_path.unlink()
                except OSError:
                    pass

        run_cubegen(cmd)

        written: Dict[int, Path] = {}
        missing: List[int] = []
        produced_for_mo: Dict[int, Path] = {}

        for mo in mos:
            mo_int = int(mo)
            suffixed_path = expected[mo_int]
            if suffixed_path.exists():
                produced_for_mo[mo_int] = suffixed_path
            elif len(mos) == 1 and tmp_path.exists():
                produced_for_mo[mo_int] = tmp_path
            else:
                missing.append(mo_int)

        if missing:
            extra_note = ""
            if len(mos) == 1:
                extra_note = f" The base single-orbital output was also checked: {tmp_path.name}"
            raise FileNotFoundError(
                "Gaussian-compatible cube generator reported success, but the expected per-orbital cube files were not found. "
                f"Missing MOs: {_format_mo_list(missing)}. Expected naming pattern: "
                f"{tmp_path.stem}_MO<N>{tmp_path.suffix}.{extra_note}"
            )

        for mo in mos:
            mo_int = int(mo)
            produced_path = produced_for_mo[mo_int]
            target_path = paths[mo_int]
            target_path.parent.mkdir(parents=True, exist_ok=True)
            if target_path.exists():
                target_path.unlink()
            produced_path.replace(target_path)
            written[mo_int] = target_path

        return written

    def _run_many_multi(items: List[tuple[int, Path, List[str]]], label: str = "  cubegen-multi-separate") -> Dict[int, Path]:
        if not items:
            return {}
        if multi_chunk <= 1:
            return _run_many(items, label=label)
        chunks = _chunk_pending(items)
        local_t0 = time.perf_counter()
        local_workers = min(jobs, len(chunks))
        results: Dict[int, Path] = {}
        if jobs == 1 or len(chunks) == 1:
            for idx, chunk in enumerate(chunks, 1):
                eprint(progress_line(label, idx - 1, len(chunks), local_t0))
                results.update(_run_one_multi(chunk))
        else:
            completed = 0
            eprint(progress_line(label, completed, len(chunks), local_t0))
            with ThreadPoolExecutor(max_workers=local_workers) as executor:
                future_to_chunk = {executor.submit(_run_one_multi, chunk): chunk[0] for chunk in chunks}
                for future in as_completed(future_to_chunk):
                    results.update(future.result())
                    completed += 1
                    eprint(progress_line(label, completed, len(chunks), local_t0))
        eprint(progress_line(label, len(chunks), len(chunks), local_t0))
        return results

    generated = _run_many_multi(pending, label="  cubegen-multi-separate") if multi_chunk > 1 else _run_many(pending, label="  cubegen")
    mo_cube_paths.update(generated)
    mo_cube_paths = _validate_or_repair_cube_headers(
        mo_cube_paths,
        generated_mos=sorted(generated),
        args=args,
        pending_builder=_build_item,
        runner=_run_many,
    )
    return mo_cube_paths
