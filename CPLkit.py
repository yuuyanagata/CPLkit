#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
transition_density_cube_combined.py

Generate both Transition Electric Dipole Moment (TEDM) and Transition Magnetic
Dipole Moment (TMDM) density cube files from a Gaussian TD-DFT output.

This script merges the TEDM workflow in tedm_density_cube_stable.py and the
TMDM workflow in tmdm_density_cube.py. It reuses a common MO-cube generation
and parsing layer, and produces TEDM and TMDM cube files in one execution.

Implemented formulas:
TEDM:
  rho_x^mu(r) = sum_j w_j * phi_i(r) * x(r) * phi_f(r)
  rho_y^mu(r) = sum_j w_j * phi_i(r) * y(r) * phi_f(r)
  rho_z^mu(r) = sum_j w_j * phi_i(r) * z(r) * phi_f(r)
  rho_total^mu(r) = (mu_x/|mu|) rho_x + (mu_y/|mu|) rho_y + (mu_z/|mu|) rho_z

TMDM:
  rho^m(r) is constructed from the orbital-gradient expression based on
  r x grad(phi), following the implementation logic in tmdm_density_cube.py:
    A(phi) = r x grad(phi)
    rho_k^m(r) = sum_j [0.5 * w_j * (phi_f * A_k(phi_i) - phi_i * A_k(phi_f))]
  rho_total^m(r) = (m_x/|m|) rho_x + (m_y/|m|) rho_y + (m_z/|m|) rho_z

Weights:
  w_j = 2 * c_j * abs(c_j)

Notes:
- Preserves Gaussian cube nvals / dataset-id line for GaussView compatibility.
- Handles natoms < 0 in cube headers correctly.
- Supports either --chk with cubegen or --mocubes_dir with pre-generated cubes.
- Uses robust cubegen lookup from the TEDM script.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import shlex
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

BOHR_PER_ANGSTROM = 1.0 / 0.52917721092
ELECTRIC_CONST_CGS = 254.17
MAGNETIC_CONST_CGS = -0.927401
HARTREE_TO_KJMOL = 2625.5
EV_TO_KJMOL = 96.4853


class Timer:
    def __init__(self) -> None:
        self.t0 = time.perf_counter()

    def elapsed(self) -> float:
        return time.perf_counter() - self.t0

    @staticmethod
    def fmt(seconds: float) -> str:
        if seconds < 60:
            return f"{seconds:5.1f}s"
        if seconds < 3600:
            m = int(seconds // 60)
            s = seconds - 60 * m
            return f"{m:d}m{s:04.1f}s"
        h = int(seconds // 3600)
        m = int((seconds - 3600 * h) // 60)
        s = seconds - 3600 * h - 60 * m
        return f"{h:d}h{m:02d}m{s:04.1f}s"


def progress_line(prefix: str, i: int, n: int, t0: float) -> str:
    now = time.perf_counter()
    done = i / max(n, 1)
    elapsed = now - t0
    eta = (elapsed / done - elapsed) if done > 0 else float("inf")
    bar_w = 26
    fill = int(bar_w * done)
    bar = "[" + "#" * fill + "-" * (bar_w - fill) + "]"
    eta_s = "  ETA " + Timer.fmt(eta) if math.isfinite(eta) else ""
    return f"{prefix} {bar} {i}/{n}  elapsed {Timer.fmt(elapsed)}{eta_s}"


def eprint(msg: str) -> None:
    print(msg, file=sys.stderr, flush=True)


@dataclass(frozen=True)
class Excitation:
    occ: int
    virt: int
    coeff: float


@dataclass
class Cube:
    comment1: str
    comment2: str
    natoms_raw: int
    natoms: int
    nvals: Optional[int]
    val_ids: Optional[List[int]]
    origin_bohr: np.ndarray
    nx: int
    ny: int
    nz: int
    vx_bohr: np.ndarray
    vy_bohr: np.ndarray
    vz_bohr: np.ndarray
    atoms_block: List[str]
    data: np.ndarray
    unit_scale_to_bohr: float


def parse_excited_state_block(log_text: str, state: int) -> List[Excitation]:
    pat_start = re.compile(rf"^\s*Excited State\s+{state}\s*:", re.MULTILINE)
    m = pat_start.search(log_text)
    if not m:
        raise ValueError(f"Could not find 'Excited State {state}:' in the log.")
    start = m.end()

    m2 = re.search(r"^\s*Excited State\s+\d+\s*:", log_text[start:], flags=re.MULTILINE)
    end = start + m2.start() if m2 else len(log_text)
    block = log_text[start:end]

    excitations: List[Excitation] = []
    line_re = re.compile(
        r"^\s*(\d+)\s*->\s*(\d+)\s+([+-]?\d+\.\d+(?:[Ee][+-]?\d+)?)",
        re.MULTILINE,
    )
    for mo_i, mo_f, coeff in line_re.findall(block):
        excitations.append(Excitation(int(mo_i), int(mo_f), float(coeff)))

    if not excitations:
        raise ValueError(
            f"Found Excited State {state} header, but no MO transitions were parsed in that block."
        )
    return excitations


def parse_tedm_vector(log_text: str, state: int) -> Tuple[float, float, float]:
    header_pat = re.compile(
        r"Ground to excited state transition electric dipole moments\s*\(Au\):",
        re.IGNORECASE,
    )
    m = header_pat.search(log_text)
    if not m:
        header_pat2 = re.compile(r"transition electric dipole moments\s*\(Au\):", re.IGNORECASE)
        m = header_pat2.search(log_text)
    if not m:
        raise ValueError("Could not find the transition electric dipole moments (Au) table in the log.")

    tail = log_text[m.end():]
    row_pat = re.compile(
        rf"^\s*{state}\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)",
        re.MULTILINE,
    )
    row = row_pat.search(tail)
    if not row:
        raise ValueError(f"Could not find TEDM vector row for state {state}.")
    return tuple(map(float, row.groups()[:3]))  # type: ignore[return-value]


def parse_tmdm_vector(log_text: str, state: int) -> Tuple[float, float, float]:
    header_pat = re.compile(
        r"Ground to excited state transition magnetic dipole moments\s*\(Au\):",
        re.IGNORECASE,
    )
    m = header_pat.search(log_text)
    if not m:
        raise ValueError("Could not find the transition magnetic dipole moments (Au) table in the log.")
    tail = log_text[m.end():]

    row_pat = re.compile(
        rf"^\s*{state}\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)",
        re.MULTILINE,
    )
    row = row_pat.search(tail)
    if not row:
        raise ValueError(f"Could not find TMDM vector row for state {state}.")
    return tuple(map(float, row.groups()[:3]))  # type: ignore[return-value]


@dataclass(frozen=True)
class CPLRow:
    state: int
    electric_x: float
    electric_y: float
    electric_z: float
    magnetic_x: float
    magnetic_y: float
    magnetic_z: float
    mu: float
    m: float
    cos_theta: float
    e_ev: float
    wavelength_nm: float
    g_value: float
    oscillator_strength: float
    s1_energy_kjmol: float
    s0_energy_kjmol: float


def _find_last_header(log_text: str, pattern: str) -> re.Match[str]:
    matches = list(re.finditer(pattern, log_text, flags=re.IGNORECASE | re.MULTILINE))
    if not matches:
        raise ValueError(f"Could not find section header matching: {pattern}")
    return matches[-1]


def _slice_section_after_header(log_text: str, header_pattern: str, end_patterns: Sequence[str]) -> str:
    header = _find_last_header(log_text, header_pattern)
    tail = log_text[header.end():]
    end_positions: List[int] = []
    for pat in end_patterns:
        m = re.search(pat, tail, flags=re.IGNORECASE | re.MULTILINE)
        if m:
            end_positions.append(m.start())
    end = min(end_positions) if end_positions else len(tail)
    return tail[:end]


def parse_tedm_table_all_states(log_text: str) -> Dict[int, Tuple[float, float, float, float]]:
    block = _slice_section_after_header(
        log_text,
        r"^\s*Ground to excited state transition electric dipole moments\s*\(Au\):",
        [
            r"^\s*Ground to excited state transition velocity dipole moments",
            r"^\s*Ground to excited state transition magnetic dipole moments",
        ],
    )
    rows: Dict[int, Tuple[float, float, float, float]] = {}
    row_re = re.compile(
        r"^\s*(\d+)\s+([+-]?\d*\.?\d+(?:[Ee][+-]?\d+)?)\s+([+-]?\d*\.?\d+(?:[Ee][+-]?\d+)?)\s+([+-]?\d*\.?\d+(?:[Ee][+-]?\d+)?)\s+([+-]?\d*\.?\d+(?:[Ee][+-]?\d+)?)\s+([+-]?\d*\.?\d+(?:[Ee][+-]?\d+)?)",
        re.MULTILINE,
    )
    for m in row_re.finditer(block):
        state = int(m.group(1))
        x, y, z, dip_s, osc = map(float, m.groups()[1:6])
        rows[state] = (x, y, z, osc)
    if not rows:
        raise ValueError("Could not parse any rows from the transition electric dipole table.")
    return rows


def parse_tmdm_table_all_states(log_text: str) -> Dict[int, Tuple[float, float, float]]:
    block = _slice_section_after_header(
        log_text,
        r"^\s*Ground to excited state transition magnetic dipole moments\s*\(Au\):",
        [
            r"^\s*Ground to excited state transition velocity quadrupole moments",
            r"^\s*Excitation energies and oscillator strengths",
        ],
    )
    rows: Dict[int, Tuple[float, float, float]] = {}
    row_re = re.compile(
        r"^\s*(\d+)\s+([+-]?\d*\.?\d+(?:[Ee][+-]?\d+)?)\s+([+-]?\d*\.?\d+(?:[Ee][+-]?\d+)?)\s+([+-]?\d*\.?\d+(?:[Ee][+-]?\d+)?)",
        re.MULTILINE,
    )
    for m in row_re.finditer(block):
        state = int(m.group(1))
        x, y, z = map(float, m.groups()[1:4])
        rows[state] = (x, y, z)
    if not rows:
        raise ValueError("Could not parse any rows from the transition magnetic dipole table.")
    return rows


def parse_excited_state_summary_all_states(log_text: str) -> Dict[int, Tuple[float, float, float]]:
    rows: Dict[int, Tuple[float, float, float]] = {}
    row_re = re.compile(
        r"^\s*Excited State\s+(\d+)\s*:\s*.*?([+-]?\d*\.?\d+(?:[Ee][+-]?\d+)?)\s+eV\s+([+-]?\d*\.?\d+(?:[Ee][+-]?\d+)?)\s+nm\s+f=([+-]?\d*\.?\d+(?:[Ee][+-]?\d+)?)",
        re.MULTILINE,
    )
    for m in row_re.finditer(log_text):
        state = int(m.group(1))
        e_ev = float(m.group(2))
        wavelength_nm = float(m.group(3))
        osc = float(m.group(4))
        rows[state] = (e_ev, wavelength_nm, osc)
    if not rows:
        raise ValueError("Could not parse any excited-state summary rows.")
    return rows


def parse_ground_state_energy_kjmol(log_text: str) -> float:
    matches = list(
        re.finditer(
            r"SCF Done:\s+E\([^)]+\)\s*=\s*([+-]?\d*\.?\d+(?:[Ee][+-]?\d+)?)",
            log_text,
            flags=re.IGNORECASE,
        )
    )
    if not matches:
        raise ValueError("Could not find any 'SCF Done' energy in the log.")
    s0_hartree = float(matches[-1].group(1))
    return s0_hartree * HARTREE_TO_KJMOL


def build_cpl_rows(log_text: str) -> List[CPLRow]:
    electric_rows = parse_tedm_table_all_states(log_text)
    magnetic_rows = parse_tmdm_table_all_states(log_text)
    excited_rows = parse_excited_state_summary_all_states(log_text)
    s0_energy_kjmol = parse_ground_state_energy_kjmol(log_text)

    states = sorted(set(electric_rows) & set(magnetic_rows) & set(excited_rows))
    if not states:
        raise ValueError("No common excited states were found across electric, magnetic, and excitation tables.")

    rows: List[CPLRow] = []
    for state in states:
        ex_au, ey_au, ez_au, osc_from_electric = electric_rows[state]
        mx_au, my_au, mz_au = magnetic_rows[state]
        e_ev, wavelength_nm, osc_from_excited = excited_rows[state]

        electric_x = ex_au * ELECTRIC_CONST_CGS
        electric_y = ey_au * ELECTRIC_CONST_CGS
        electric_z = ez_au * ELECTRIC_CONST_CGS
        magnetic_x = mx_au * MAGNETIC_CONST_CGS
        magnetic_y = my_au * MAGNETIC_CONST_CGS
        magnetic_z = mz_au * MAGNETIC_CONST_CGS

        mu = math.sqrt(electric_x * electric_x + electric_y * electric_y + electric_z * electric_z)
        m = math.sqrt(magnetic_x * magnetic_x + magnetic_y * magnetic_y + magnetic_z * magnetic_z)
        dot = electric_x * magnetic_x + electric_y * magnetic_y + electric_z * magnetic_z
        cos_theta = dot / (mu * m) if mu > 0.0 and m > 0.0 else float("nan")
        g_value = (4.0 * mu * m * cos_theta / (mu * mu + m * m)) if mu > 0.0 and m > 0.0 else float("nan")
        oscillator_strength = osc_from_excited
        if math.isfinite(osc_from_electric) and abs(osc_from_electric - osc_from_excited) > 1e-6:
            oscillator_strength = osc_from_excited
        s1_energy_kjmol = s0_energy_kjmol + e_ev * EV_TO_KJMOL

        rows.append(
            CPLRow(
                state=state,
                electric_x=electric_x,
                electric_y=electric_y,
                electric_z=electric_z,
                magnetic_x=magnetic_x,
                magnetic_y=magnetic_y,
                magnetic_z=magnetic_z,
                mu=mu,
                m=m,
                cos_theta=cos_theta,
                e_ev=e_ev,
                wavelength_nm=wavelength_nm,
                g_value=g_value,
                oscillator_strength=oscillator_strength,
                s1_energy_kjmol=s1_energy_kjmol,
                s0_energy_kjmol=s0_energy_kjmol,
            )
        )
    return rows


def write_cpl_csv(path: Path, source_log: Path, rows: Sequence[CPLRow]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    header = [
        "file",
        "state",
        "electric X",
        "electric Y",
        "electric Z",
        "magnetic X",
        "magnetic Y",
        "magnetic Z",
        "mu",
        "m",
        "cos theta",
        "eV",
        "nm",
        "g value",
        "Oscillator Strength",
        "S1-Energy",
        "S0-Energy",
    ]
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        for row in rows:
            writer.writerow(
                [
                    source_log.stem,
                    row.state,
                    f"{row.electric_x:.18g}",
                    f"{row.electric_y:.18g}",
                    f"{row.electric_z:.18g}",
                    f"{row.magnetic_x:.18g}",
                    f"{row.magnetic_y:.18g}",
                    f"{row.magnetic_z:.18g}",
                    f"{row.mu:.18g}",
                    f"{row.m:.18g}",
                    f"{row.cos_theta:.18g}",
                    f"{row.e_ev:.18g}",
                    f"{row.wavelength_nm:.18g}",
                    f"{row.g_value:.18g}",
                    f"{row.oscillator_strength:.18g}",
                    f"{row.s1_energy_kjmol:.18g}",
                    f"{row.s0_energy_kjmol:.18g}",
                ]
            )


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


def read_cube(path: Path, dtype: np.dtype) -> Cube:
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

    unit_scale_to_bohr = BOHR_PER_ANGSTROM if natoms_raw < 0 else 1.0

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
    arr = np.fromstring(data_text, sep=" ", dtype=float)
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


def orbital_gradients(phi: np.ndarray, cube: Cube) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    dx, dy, dz = cube_step_sizes_bohr(cube)
    dphidx, dphidy, dphidz = np.gradient(phi, dx, dy, dz, edge_order=2)
    return dphidx, dphidy, dphidz


def precompute_r_cross_grad(
    cube: Cube,
    grad: Tuple[np.ndarray, np.ndarray, np.ndarray],
    coords_mode: str,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    dphidx, dphidy, dphidz = grad

    if coords_mode == "aligned":
        x, y, z = coordinate_1d_aligned_bohr(cube)
        ax = (y[None, :, None] * dphidz) - (z[None, None, :] * dphidy)
        ay = (z[None, None, :] * dphidx) - (x[:, None, None] * dphidz)
        az = (x[:, None, None] * dphidy) - (y[None, :, None] * dphidx)
        return ax, ay, az

    x, y, z = coordinate_3d_general_bohr(cube)
    ax = y * dphidz - z * dphidy
    ay = z * dphidx - x * dphidz
    az = x * dphidy - y * dphidx
    return ax, ay, az


def excitation_weight(coeff: float) -> float:
    return 2.0 * coeff * abs(coeff)


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
) -> List[str]:
    cmd = [cubegen_exe, str(npts), f"MO={mo}", str(chk_path), str(cube_path)]
    cmd.extend(grid_args)
    return cmd


def run_cubegen(cmd: Sequence[str]) -> None:
    res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if res.returncode != 0:
        msg = [
            f"cubegen failed with exit status {res.returncode}.",
            "Command:",
            "  " + subprocess.list2cmdline(list(cmd)),
        ]
        if res.stdout and res.stdout.strip():
            msg.append("stdout:\n" + res.stdout.rstrip())
        if res.stderr and res.stderr.strip():
            msg.append("stderr:\n" + res.stderr.rstrip())
        if os.name == "nt" and res.returncode == 127:
            msg.append(
                "On Windows, exit status 127 often means that cubegen or a dependent Gaussian runtime DLL was not found."
            )
        raise RuntimeError("\n".join(msg))


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
        return mo_cube_paths

    if args.chk is None:
        raise ValueError("Provide either --chk or --mocubes_dir.")

    cubegen_exe = resolve_cubegen(args.cubegen)
    grid_args = shlex.split(args.cubegen_grid, posix=False)
    mocdir = outdir / f"{baseprefix}_mocubes"
    mocdir.mkdir(parents=True, exist_ok=True)

    eprint("[3/8] Generating MO cubes with cubegen ...")
    eprint(f"[Info] cubegen = {cubegen_exe}")
    gen_t0 = time.perf_counter()

    for idx, mo in enumerate(needed_mos, 1):
        eprint(progress_line("  cubegen", idx - 1, len(needed_mos), gen_t0))
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
        )
        run_cubegen(cmd)
        if not cube_path.exists():
            raise FileNotFoundError(f"cubegen reported success but output cube was not created: {cube_path}")
        mo_cube_paths[mo] = cube_path

    eprint(progress_line("  cubegen", len(needed_mos), len(needed_mos), gen_t0))
    return mo_cube_paths


def main() -> None:
    ap = argparse.ArgumentParser(
        description=(
            "Generate TEDM and TMDM density cube files from Gaussian TDDFT output and "
            "optionally export CPL g values to CSV."
        )
    )
    ap.add_argument("--log", required=True, type=Path, help="Gaussian TDDFT output (*.log/*.out)")
    ap.add_argument("--state", type=int, default=None, help="Excited state index (1-based) for cube generation")
    ap.add_argument("--chk", type=Path, default=None, help="Checkpoint (chk/fchk) for cubegen")
    ap.add_argument("--mocubes_dir", type=Path, default=None, help="Directory with mo<MO>.cube files")
    ap.add_argument("--cubegen", default="cubegen", help="Path to cubegen executable")
    ap.add_argument("--cubegen_npts", type=int, default=0, help="First cubegen arg")
    ap.add_argument("--cubegen_grid", default="-3 h", help='Remaining cubegen args, for example "-3 h"')
    ap.add_argument("--outdir", type=Path, default=Path("."), help="Output directory")
    ap.add_argument("--outprefix", default=None, help="Output prefix. Default is S<state>")
    ap.add_argument("--overwrite_mo_cubes", action="store_true", help="Force regeneration of MO cubes")
    ap.add_argument("--keep_components", action="store_true", help="Also write x, y, z component cubes")
    ap.add_argument(
        "--dtype",
        choices=["float32", "float64"],
        default="float32",
        help="Internal dtype",
    )
    ap.add_argument(
        "--coords",
        choices=["auto", "aligned", "general"],
        default="auto",
        help="Coordinate handling mode",
    )
    ap.add_argument(
        "--cpl_only",
        action="store_true",
        help="Only export the CPL CSV and skip cube generation",
    )
    ap.add_argument(
        "--no_cpl_csv",
        action="store_true",
        help="Disable CPL CSV export",
    )
    ap.add_argument(
        "--cpl_csv_path",
        type=Path,
        default=None,
        help="Output path for the CPL CSV. Default is <logstem>-CPL.csv in --outdir",
    )
    args = ap.parse_args()

    total_timer = Timer()
    eprint(f"[Start] {args.log.name}  state={args.state}")

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
    mux, muy, muz = parse_tedm_vector(log_text, args.state)
    mx, my, mz = parse_tmdm_vector(log_text, args.state)
    mu_norm = math.sqrt(mux * mux + muy * muy + muz * muz)
    m_norm = math.sqrt(mx * mx + my * my + mz * mz)
    if mu_norm == 0.0:
        raise ValueError(f"TEDM vector norm is zero for state {args.state}.")
    if m_norm == 0.0:
        raise ValueError(f"TMDM vector norm is zero for state {args.state}.")
    eprint(
        f"[2/8] Parsed target state in {Timer.fmt(time.perf_counter() - t0)}: "
        f"excitations={len(excitations)} |mu|={mu_norm:.6f} Au |m|={m_norm:.6f} Au"
    )

    outprefix = args.outprefix or f"S{args.state}"
    logstem = args.log.stem
    baseprefix = f"{logstem}_{outprefix}"
    needed_mos = sorted({ex.occ for ex in excitations} | {ex.virt for ex in excitations})
    eprint(f"[Info] Needed MOs: {needed_mos}")
    mo_cube_paths = ensure_mo_cubes(needed_mos, args, baseprefix, outdir)

    t0 = time.perf_counter()
    ref_cube = read_cube(mo_cube_paths[needed_mos[0]], dtype=dtype)
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

    if coords_mode == "aligned":
        x1d, y1d, z1d = coordinate_1d_aligned_bohr(ref_cube)
        x = x1d[:, None, None].astype(dtype, copy=False)
        y = y1d[None, :, None].astype(dtype, copy=False)
        z = z1d[None, None, :].astype(dtype, copy=False)
    else:
        x, y, z = coordinate_3d_general_bohr(ref_cube)
        x = x.astype(dtype, copy=False)
        y = y.astype(dtype, copy=False)
        z = z.astype(dtype, copy=False)

    phi_cache: Dict[int, np.ndarray] = {}
    a_cache: Dict[int, Tuple[np.ndarray, np.ndarray, np.ndarray]] = {}

    def get_phi(mo: int) -> np.ndarray:
        if mo not in phi_cache:
            c = read_cube(mo_cube_paths[mo], dtype=dtype)
            if (c.nx, c.ny, c.nz) != (ref_cube.nx, ref_cube.ny, ref_cube.nz):
                raise ValueError(f"Grid mismatch for mo{mo}")
            if not np.isclose(c.unit_scale_to_bohr, ref_cube.unit_scale_to_bohr):
                raise ValueError(f"Unit conversion mismatch for mo{mo}. Regenerate cubes consistently.")
            if c.nvals != ref_cube.nvals:
                raise ValueError(f"nvals mismatch for mo{mo}. Regenerate cubes consistently.")
            phi_cache[mo] = c.data
        return phi_cache[mo]

    def get_A(mo: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        if mo in a_cache:
            return a_cache[mo]
        phi = get_phi(mo)
        grad = orbital_gradients(phi, ref_cube)
        ax, ay, az = precompute_r_cross_grad(ref_cube, grad, coords_mode=coords_mode)
        a_cache[mo] = (
            np.asarray(ax, dtype=dtype),
            np.asarray(ay, dtype=dtype),
            np.asarray(az, dtype=dtype),
        )
        return a_cache[mo]

    rho_tedm_x = np.zeros((ref_cube.nx, ref_cube.ny, ref_cube.nz), dtype=dtype)
    rho_tedm_y = np.zeros_like(rho_tedm_x)
    rho_tedm_z = np.zeros_like(rho_tedm_x)

    rho_tmdm_x = np.zeros_like(rho_tedm_x)
    rho_tmdm_y = np.zeros_like(rho_tedm_x)
    rho_tmdm_z = np.zeros_like(rho_tedm_x)

    eprint("[5/8] Accumulating excitations for TEDM and TMDM ...")
    ex_t0 = time.perf_counter()
    n_ex = len(excitations)
    for i, ex in enumerate(excitations, 1):
        if i == 1 or i % 5 == 0 or i == n_ex:
            eprint(progress_line("  excit", i - 1, n_ex, ex_t0))

        w = excitation_weight(ex.coeff)
        phi_occ = get_phi(ex.occ)
        phi_virt = get_phi(ex.virt)
        tau = phi_occ * phi_virt

        rho_tedm_x += w * (x * tau)
        rho_tedm_y += w * (y * tau)
        rho_tedm_z += w * (z * tau)

        s = 0.5 * w
        a_occ = get_A(ex.occ)
        a_virt = get_A(ex.virt)
        rho_tmdm_x += s * (phi_virt * a_occ[0] - phi_occ * a_virt[0])
        rho_tmdm_y += s * (phi_virt * a_occ[1] - phi_occ * a_virt[1])
        rho_tmdm_z += s * (phi_virt * a_occ[2] - phi_occ * a_virt[2])

    eprint(progress_line("  excit", n_ex, n_ex, ex_t0))

    eprint("[6/8] Projecting onto transition vectors ...")
    muxn, muyn, muzn = mux / mu_norm, muy / mu_norm, muz / mu_norm
    mxn, myn, mzn = mx / m_norm, my / m_norm, mz / m_norm

    rho_tedm_total = muxn * rho_tedm_x + muyn * rho_tedm_y + muzn * rho_tedm_z
    rho_tmdm_total = mxn * rho_tmdm_x + myn * rho_tmdm_y + mzn * rho_tmdm_z

    eprint("[7/8] Writing cube files ...")
    tedm_comment1 = f"TEDM density cube from {args.log.name} (state {args.state})"
    tedm_comment2 = f"Kubo JPCL 2021 SI eqs. S11-S13; projected onto mu (Au). dtype={args.dtype}"
    tmdm_comment1 = f"TMDM density cube from {args.log.name} (state {args.state})"
    tmdm_comment2 = f"Kubo JPCL 2021 SI eqs. S14-S20; projected onto m (Au). dtype={args.dtype}"

    tedm_total_path = outdir / f"{baseprefix}_TEDM_total.cube"
    tmdm_total_path = outdir / f"{baseprefix}_TMDM_total.cube"

    write_cube(tedm_total_path, ref_cube, rho_tedm_total, tedm_comment1, tedm_comment2)
    write_cube(tmdm_total_path, ref_cube, rho_tmdm_total, tmdm_comment1, tmdm_comment2)

    if args.keep_components:
        write_cube(outdir / f"{baseprefix}_TEDM_x.cube", ref_cube, rho_tedm_x, tedm_comment1, tedm_comment2)
        write_cube(outdir / f"{baseprefix}_TEDM_y.cube", ref_cube, rho_tedm_y, tedm_comment1, tedm_comment2)
        write_cube(outdir / f"{baseprefix}_TEDM_z.cube", ref_cube, rho_tedm_z, tedm_comment1, tedm_comment2)

        write_cube(outdir / f"{baseprefix}_TMDM_x.cube", ref_cube, rho_tmdm_x, tmdm_comment1, tmdm_comment2)
        write_cube(outdir / f"{baseprefix}_TMDM_y.cube", ref_cube, rho_tmdm_y, tmdm_comment1, tmdm_comment2)
        write_cube(outdir / f"{baseprefix}_TMDM_z.cube", ref_cube, rho_tmdm_z, tmdm_comment1, tmdm_comment2)

    eprint("[8/8] Done")
    eprint(f"  Total elapsed: {Timer.fmt(total_timer.elapsed())}")
    eprint(f"  TEDM output: {tedm_total_path}")
    eprint(f"  TMDM output: {tmdm_total_path}")
    if not args.no_cpl_csv:
        eprint(f"  CPL CSV output: {cpl_csv_path}")


if __name__ == "__main__":
    main()
