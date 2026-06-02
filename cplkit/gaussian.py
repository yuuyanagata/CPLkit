"""Gaussian TD-DFT output parsers and CPL summary calculations.

This module is part of an independently written Python implementation. It parses
Gaussian table labels literally, including labels that contain the phrase
"transition electric dipole moments", while the manuscript prose uses
the abbreviation "EDTM" for electric dipole transition moment.
"""

from __future__ import annotations

import csv
import math
import re
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from .constants import (
    ELECTRIC_CONST_CGS,
    EV_TO_KJMOL,
    HARTREE_TO_KJMOL,
    MAGNETIC_CONST_CGS,
    RADIATIVE_RATE_CONST_NS,
)
from .models import CoefficientDiagnostics, CPLRow, Excitation


def parse_excited_state_block(log_text: str, state: int) -> List[Excitation]:
    """Parse the configuration coefficients printed under one Gaussian excited-state block.

    Gaussian normally prints only configurations above an internal threshold. The returned
    list must therefore be interpreted together with :func:`diagnose_coefficients` when
    density reconstruction is intended.
    """

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


def parse_route_section(log_text: str) -> str:
    """Return the first Gaussian route section as a single lowercase string."""

    lines = log_text.splitlines()
    route: List[str] = []
    in_route = False
    for line in lines:
        stripped = line.strip()
        if not in_route and stripped.startswith("#"):
            in_route = True
        if in_route:
            if not stripped and route:
                break
            route.append(stripped)
    return " ".join(route).lower()


def detect_iop_940_level(log_text: str) -> Optional[int]:
    """Detect whether IOp(9/40=n) appears in the Gaussian route section."""

    route = parse_route_section(log_text)
    m = re.search(r"iop\s*\(\s*9\s*/\s*40\s*=\s*(\d+)\s*\)", route, flags=re.IGNORECASE)
    if not m:
        return None
    return int(m.group(1))


def diagnose_coefficients(log_text: str, state: int, excitations: Sequence[Excitation]) -> CoefficientDiagnostics:
    """Assess whether the printed CI coefficients are likely sufficient for density reconstruction.

    For closed-shell singlet TD-DFT output in Gaussian, the printed coefficients commonly satisfy
    approximately 2 * sum(c_i^2) = 1 when the dominant configurations have been retained. A much
    smaller value is a practical warning that the default Gaussian print threshold has removed
    non-negligible configurations. This diagnostic does not replace a convergence study.
    """

    two_sum_c_squared = 2.0 * sum(ex.coeff * ex.coeff for ex in excitations)
    iop_level = detect_iop_940_level(log_text)
    has_recommended_iop = iop_level is not None and iop_level >= 4

    messages: List[str] = []
    if not has_recommended_iop:
        messages.append(
            "Gaussian route section does not contain IOp(9/40=4) or a stricter equivalent; "
            "the printed TD-DFT configurations may be threshold-truncated."
        )
    if two_sum_c_squared < 0.98:
        messages.append(
            f"2*sum(c^2) = {two_sum_c_squared:.5f}, which suggests incomplete printed configurations."
        )
    if not messages:
        messages.append("No coefficient-printing warning was detected by the built-in heuristic.")

    return CoefficientDiagnostics(
        state=state,
        n_configurations=len(excitations),
        two_sum_c_squared=two_sum_c_squared,
        iop_940_level=iop_level,
        has_recommended_iop=has_recommended_iop,
        warning=" ".join(messages),
    )


def parse_edtm_vector(log_text: str, state: int) -> Tuple[float, float, float]:
    header_pat = re.compile(
        r"Ground to excited state transition electric dipole moments\s*\(Au\):",
        re.IGNORECASE,
    )
    m = header_pat.search(log_text)
    if not m:
        header_pat2 = re.compile(r"transition electric dipole moments\s*\(Au\):", re.IGNORECASE)
        m = header_pat2.search(log_text)
    if not m:
        raise ValueError("Could not find the Gaussian EDTM table in the log.")

    tail = log_text[m.end():]
    row_pat = re.compile(
        rf"^\s*{state}\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)",
        re.MULTILINE,
    )
    row = row_pat.search(tail)
    if not row:
        raise ValueError(f"Could not find EDTM vector row for state {state}.")
    return tuple(map(float, row.groups()[:3]))  # type: ignore[return-value]


def parse_mdtm_vector(log_text: str, state: int) -> Tuple[float, float, float]:
    header_pat = re.compile(
        r"Ground to excited state transition magnetic dipole moments\s*\(Au\):",
        re.IGNORECASE,
    )
    m = header_pat.search(log_text)
    if not m:
        raise ValueError("Could not find the Gaussian MDTM table in the log.")
    tail = log_text[m.end():]

    row_pat = re.compile(
        rf"^\s*{state}\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)",
        re.MULTILINE,
    )
    row = row_pat.search(tail)
    if not row:
        raise ValueError(f"Could not find MDTM vector row for state {state}.")
    return tuple(map(float, row.groups()[:3]))  # type: ignore[return-value]


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


def parse_edtm_table_all_states(log_text: str) -> Dict[int, Tuple[float, float, float, float]]:
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
        x, y, z, _dip_s, osc = map(float, m.groups()[1:6])
        rows[state] = (x, y, z, osc)
    if not rows:
        raise ValueError("Could not parse any rows from the Gaussian EDTM table.")
    return rows


def parse_mdtm_table_all_states(log_text: str) -> Dict[int, Tuple[float, float, float]]:
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
        raise ValueError("Could not parse any rows from the Gaussian MDTM table.")
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
    """Build state-wise CPL summary rows from Gaussian transition-moment tables.

    The calculation uses the EDTM and MDTM values printed by Gaussian.  The vectors are converted to
    the CGS-based units used for the CPL summary workflow.  For each excited
    state, the function evaluates |mu|, |m|, cos(theta), and

        g_CPL = 4 |mu| |m| cos(theta) / (|mu|^2 + |m|^2).

    This routine does not use cube integration.  It summarizes Gaussian tables
    directly and is therefore separate from response-density validation.
    """

    electric_rows = parse_edtm_table_all_states(log_text)
    magnetic_rows = parse_mdtm_table_all_states(log_text)
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
        radiative_rate_constant_ns = (
            RADIATIVE_RATE_CONST_NS * oscillator_strength / (wavelength_nm * wavelength_nm)
        ) if wavelength_nm > 0.0 else float("nan")
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
                radiative_rate_constant_ns=radiative_rate_constant_ns,
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
        "Radiative Rate Constant (ns^-1)",
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
                    f"{row.radiative_rate_constant_ns:.18g}",
                ]
            )
