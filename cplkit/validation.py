"""Validation helpers for cube-integrated moment-density components."""

from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np

from .cubeio import cube_voxel_volume_bohr3
from .models import Cube, MomentValidation


def compare_component_integrals(
    cube: Cube,
    quantity: str,
    component_grids: Sequence[np.ndarray],
    gaussian_components_au: Sequence[float],
) -> list[MomentValidation]:
    """Compare integrated density components with Gaussian transition-moment components."""

    vol = cube_voxel_volume_bohr3(cube)
    rows: list[MomentValidation] = []
    for comp, grid, target in zip(("x", "y", "z"), component_grids, gaussian_components_au):
        integrated = float(np.sum(grid, dtype=np.float64) * vol)
        abs_err = abs(integrated - target)
        rel_err = abs_err / abs(target) if abs(target) > 1e-14 else (0.0 if abs_err < 1e-14 else math.inf)
        rows.append(
            MomentValidation(
                quantity=quantity,
                component=comp,
                integrated_au=integrated,
                gaussian_au=float(target),
                absolute_error=abs_err,
                relative_error=rel_err,
            )
        )
    return rows


def write_validation_csv(path: Path, rows: Iterable[MomentValidation]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["quantity", "component", "integrated_au", "gaussian_au", "absolute_error", "relative_error"])
        for row in rows:
            writer.writerow(
                [
                    row.quantity,
                    row.component,
                    f"{row.integrated_au:.18g}",
                    f"{row.gaussian_au:.18g}",
                    f"{row.absolute_error:.18g}",
                    f"{row.relative_error:.18g}",
                ]
            )
