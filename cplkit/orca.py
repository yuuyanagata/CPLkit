"""ORCA TD-DFT/TDA spectrum parsers and CPL summary calculations.

The implementation reads the length-gauge absorption and CD tables printed by
ORCA.  ORCA reports the rotatory strength R alongside the magnetic transition
dipole components.  R is used directly for the signed scalar product in the
dissymmetry-factor expression, avoiding loss of precision from rounded vector
components.
"""

from __future__ import annotations

import math
import re
from typing import Dict, List, Sequence, Tuple

from .constants import (
    ELECTRIC_CONST_CGS,
    EV_TO_KJMOL,
    HARTREE_TO_KJMOL,
    MAGNETIC_CONST_CGS,
    RADIATIVE_RATE_CONST_NS,
)
from .models import CPLRow


# ORCA's printed length-gauge CD table satisfies
# R(1e-40 cgs) = mu_cgs dot (2 * m_cgs) for its tabulated atomic-unit
# transition moments.  The factor of two is therefore required to place ORCA's
# magnetic components in the same scaled-CGS convention as the printed R.
ORCA_MAGNETIC_CONST_CGS = 2.0 * MAGNETIC_CONST_CGS

_FLOAT_PATTERN = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?"
_FLOAT_RE = re.compile(_FLOAT_PATTERN)
_TARGET_STATE_RE = re.compile(r"->\s*(\d+)\s*-", re.IGNORECASE)


def _as_float(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def _find_last_header(output_text: str, pattern: str) -> re.Match[str]:
    matches = list(re.finditer(pattern, output_text, flags=re.IGNORECASE | re.MULTILINE))
    if not matches:
        raise ValueError(f"Could not find ORCA section header matching: {pattern}")
    return matches[-1]


def _slice_last_section(
    output_text: str,
    header_pattern: str,
    end_patterns: Sequence[str],
) -> str:
    header = _find_last_header(output_text, header_pattern)
    tail = output_text[header.end():]
    end_positions: List[int] = []
    for pattern in end_patterns:
        match = re.search(pattern, tail, flags=re.IGNORECASE | re.MULTILINE)
        if match:
            end_positions.append(match.start())
    end = min(end_positions) if end_positions else len(tail)
    return tail[:end]


def _parse_transition_rows(block: str, number_count: int) -> Dict[int, Tuple[float, ...]]:
    rows: Dict[int, Tuple[float, ...]] = {}
    for line in block.splitlines():
        state_match = _TARGET_STATE_RE.search(line)
        if state_match is None:
            continue
        values = [_as_float(value) for value in _FLOAT_RE.findall(line)]
        if len(values) < number_count:
            continue
        rows[int(state_match.group(1))] = tuple(values[-number_count:])
    return rows


def parse_absorption_spectrum_all_states(
    output_text: str,
) -> Dict[int, Tuple[float, float, float, float, float, float]]:
    """Parse ORCA's last length-gauge absorption table.

    Each value tuple contains excitation energy (eV), wavelength (nm),
    oscillator strength, and the x/y/z electric transition dipole in a.u.
    """

    block = _slice_last_section(
        output_text,
        r"^\s*ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS\s*$",
        [
            r"^\s*ABSORPTION SPECTRUM VIA TRANSITION VELOCITY DIPOLE MOMENTS\s*$",
            r"^\s*CD SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS\s*$",
        ],
    )
    parsed = _parse_transition_rows(block, number_count=8)
    rows: Dict[int, Tuple[float, float, float, float, float, float]] = {}
    for state, values in parsed.items():
        e_ev, _energy_cm, wavelength_nm, oscillator, _d2, dx, dy, dz = values
        rows[state] = (e_ev, wavelength_nm, oscillator, dx, dy, dz)
    if not rows:
        raise ValueError("Could not parse any rows from ORCA's length-gauge absorption spectrum.")
    return rows


def parse_cd_spectrum_all_states(
    output_text: str,
) -> Dict[int, Tuple[float, float, float, float]]:
    """Parse ORCA's last length-gauge CD table.

    Each value tuple contains R in ORCA's printed ``1e-40 cgs`` unit followed
    by the x/y/z magnetic transition dipole in a.u.
    """

    block = _slice_last_section(
        output_text,
        r"^\s*CD SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS\s*$",
        [
            r"^\s*CD SPECTRUM VIA TRANSITION VELOCITY DIPOLE MOMENTS\s*$",
            r"^\s*CIS/TD-DFT TOTAL ENERGY\s*$",
            r"^\s*ORCA-CIS/TD-DFT FINISHED\s*$",
            r"^\s*Maximum memory used throughout the entire CIS-calculation",
        ],
    )
    parsed = _parse_transition_rows(block, number_count=7)
    rows: Dict[int, Tuple[float, float, float, float]] = {}
    for state, values in parsed.items():
        _e_ev, _energy_cm, _wavelength_nm, rotatory_strength, mx, my, mz = values
        rows[state] = (rotatory_strength, mx, my, mz)
    if not rows:
        raise ValueError("Could not parse any rows from ORCA's length-gauge CD spectrum.")
    return rows


def parse_ground_state_energy_kjmol(output_text: str) -> float:
    """Return the last ORCA ``E(SCF)`` value in kJ mol-1."""

    matches = list(
        re.finditer(
            rf"^\s*E\(SCF\)\s*=\s*({_FLOAT_PATTERN})\s+Eh\b",
            output_text,
            flags=re.IGNORECASE | re.MULTILINE,
        )
    )
    if not matches:
        raise ValueError("Could not find an ORCA 'E(SCF)' energy in the output.")
    return _as_float(matches[-1].group(1)) * HARTREE_TO_KJMOL


def build_cpl_rows(output_text: str) -> List[CPLRow]:
    """Build state-wise CPL rows from ORCA TD-DFT or TDA spectrum tables.

    The dimensionless dissymmetry factor is evaluated in the same scaled-CGS
    convention used by CPLkit's Gaussian summary:

        g_CPL = 4 R / (|mu|^2 + |m|^2)

    where ORCA's printed length-gauge R is the signed ``mu dot m`` scalar in
    ``1e-40 cgs`` after applying ORCA's magnetic-moment convention.
    """

    electric_rows = parse_absorption_spectrum_all_states(output_text)
    magnetic_rows = parse_cd_spectrum_all_states(output_text)
    try:
        s0_energy_kjmol = parse_ground_state_energy_kjmol(output_text)
    except ValueError:
        # S0/S1 total energies are not part of the exported CPL CSV.  Keep CPL
        # analysis available for otherwise complete spectrum-only outputs.
        s0_energy_kjmol = float("nan")

    states = sorted(set(electric_rows) & set(magnetic_rows))
    if not states:
        raise ValueError("No common excited states were found in ORCA's absorption and CD tables.")

    rows: List[CPLRow] = []
    for state in states:
        e_ev, wavelength_nm, oscillator_strength, ex_au, ey_au, ez_au = electric_rows[state]
        rotatory_strength, mx_au, my_au, mz_au = magnetic_rows[state]

        electric_x = ex_au * ELECTRIC_CONST_CGS
        electric_y = ey_au * ELECTRIC_CONST_CGS
        electric_z = ez_au * ELECTRIC_CONST_CGS
        magnetic_x = mx_au * ORCA_MAGNETIC_CONST_CGS
        magnetic_y = my_au * ORCA_MAGNETIC_CONST_CGS
        magnetic_z = mz_au * ORCA_MAGNETIC_CONST_CGS

        mu = math.sqrt(electric_x * electric_x + electric_y * electric_y + electric_z * electric_z)
        m = math.sqrt(magnetic_x * magnetic_x + magnetic_y * magnetic_y + magnetic_z * magnetic_z)
        if mu > 0.0 and m > 0.0:
            vector_dot = (
                electric_x * magnetic_x
                + electric_y * magnetic_y
                + electric_z * magnetic_z
            )
            cos_theta = vector_dot / (mu * m)
            g_value = 4.0 * rotatory_strength / (mu * mu + m * m)
        else:
            cos_theta = float("nan")
            g_value = float("nan")

        radiative_rate_constant_ns = (
            RADIATIVE_RATE_CONST_NS * oscillator_strength / (wavelength_nm * wavelength_nm)
            if wavelength_nm > 0.0
            else float("nan")
        )
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
