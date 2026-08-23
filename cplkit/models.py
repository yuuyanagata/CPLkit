"""Typed data containers for CPL summaries, Gaussian cube files, and validation."""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Tuple

import numpy as np


@dataclass(frozen=True)
class Excitation:
    """One Gaussian TD-DFT configuration line of the form occupied -> virtual coefficient."""

    occ: int
    virt: int
    coeff: float


@dataclass
class Cube:
    """Gaussian cube file content with numerical grid coordinates expressed internally in Bohr."""

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


@dataclass(frozen=True)
class CPLRow:
    """One row in the state-specific CPL summary CSV."""

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
    radiative_rate_constant_ns: float
    s1_energy_kjmol: float
    s0_energy_kjmol: float


@dataclass(frozen=True)
class CoefficientDiagnostics:
    """Diagnostics for the TD-DFT configuration coefficients printed by Gaussian."""

    state: int
    n_configurations: int
    two_sum_c_squared: float
    iop_940_level: Optional[int]
    has_recommended_iop: bool
    warning: str


@dataclass(frozen=True)
class MomentValidation:
    """Numerical comparison between cube-integrated transition moment densities and Gaussian tables."""

    quantity: str
    component: str
    integrated_au: float
    gaussian_au: float
    absolute_error: float
    relative_error: float
