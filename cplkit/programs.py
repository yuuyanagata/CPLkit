"""Electronic-structure program detection for CPL summary parsing."""

from __future__ import annotations

import re


SUPPORTED_PROGRAMS = ("gaussian", "orca")


def detect_output_program(output_text: str) -> str:
    """Detect whether *output_text* was written by Gaussian or ORCA.

    The transition-spectrum headers are checked before the more general program
    banners because output files are sometimes concatenated with scheduler logs.
    """

    if re.search(
        r"^\s*(?:ABSORPTION|CD) SPECTRUM VIA TRANSITION (?:ELECTRIC|VELOCITY) DIPOLE MOMENTS\s*$",
        output_text,
        flags=re.IGNORECASE | re.MULTILINE,
    ):
        return "orca"
    if re.search(
        r"Ground to excited state transition (?:electric|magnetic) dipole moments",
        output_text,
        flags=re.IGNORECASE,
    ):
        return "gaussian"

    if re.search(r"\bO\s+R\s+C\s+A\b|\bORCA\s+(?:VERSION|TERMINATED)\b", output_text, re.IGNORECASE):
        return "orca"
    if re.search(r"Gaussian, Inc\.|Entering Gaussian System|SCF Done:", output_text, re.IGNORECASE):
        return "gaussian"

    raise ValueError(
        "Could not identify the output as Gaussian or ORCA. "
        "Specify --program gaussian or --program orca explicitly."
    )
