import csv
import math

from cplkit.cli import main
from cplkit.constants import ELECTRIC_CONST_CGS, MAGNETIC_CONST_CGS
from cplkit.orca import (
    build_cpl_rows,
    parse_absorption_spectrum_all_states,
    parse_cd_spectrum_all_states,
)
from cplkit.programs import detect_output_program


ORCA_OUTPUT = """
                              O   R   C   A

                     ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS
----------------------------------------------------------------------------------------------------
     Transition      Energy     Energy  Wavelength fosc(D2)      D2        DX        DY        DZ
                      (eV)      (cm-1)    (nm)                 (au**2)    (au)      (au)      (au)
----------------------------------------------------------------------------------------------------
  0-1 Ag ->  1-1 B3u  3.811235   30739.7   325.3   0.339291983   3.63371   1.78104  -0.17825  -0.65562
  0-1 Ag ->  2-1 B2u  4.042790   32607.3   306.7   1.154479791  11.65593   2.35419   2.37443   0.68977

                     ABSORPTION SPECTRUM VIA TRANSITION VELOCITY DIPOLE MOMENTS
----------------------------------------------------------------------------------------------------

                     CD SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS
------------------------------------------------------------------------------------------
     Transition      Energy     Energy  Wavelength    R        MX        MY        MZ
                      (eV)      (cm-1)    (nm)   (1e40*cgs)   (au)      (au)      (au)
------------------------------------------------------------------------------------------
  0-1 Ag ->  1-1 B3u  3.811235   30739.7   325.3  -420.54720   0.25960  -0.03543  -0.64573
  0-1 Ag ->  2-1 B2u  4.042790   32607.3   306.7   744.90418  -0.23136   0.04506  -1.65613

                     CD SPECTRUM VIA TRANSITION VELOCITY DIPOLE MOMENTS
------------------------------------------------------------------------------------------

    E(SCF)  =  -1800.679203610 Eh
"""


def test_detects_orca_output():
    assert detect_output_program(ORCA_OUTPUT) == "orca"


def test_detects_gaussian_output():
    gaussian_output = "Ground to excited state transition electric dipole moments (Au):"
    assert detect_output_program(gaussian_output) == "gaussian"


def test_parses_orca_length_gauge_tables_with_spaced_symmetry_labels():
    absorption = parse_absorption_spectrum_all_states(ORCA_OUTPUT)
    cd = parse_cd_spectrum_all_states(ORCA_OUTPUT)

    assert absorption[1] == (3.811235, 325.3, 0.339291983, 1.78104, -0.17825, -0.65562)
    assert cd[2] == (744.90418, -0.23136, 0.04506, -1.65613)


def test_builds_orca_cpl_g_value_from_printed_rotatory_strength():
    row = build_cpl_rows(ORCA_OUTPUT)[0]

    assert row.state == 1
    assert math.isclose(row.electric_x, 1.78104 * ELECTRIC_CONST_CGS)
    assert math.isclose(row.magnetic_x, 0.25960 * 2.0 * MAGNETIC_CONST_CGS)

    expected_g = 4.0 * -420.54720 / (row.mu * row.mu + row.m * row.m)
    assert math.isclose(row.g_value, expected_g, rel_tol=1e-12)
    assert -1.0 <= row.cos_theta <= 1.0

    # The scaled vectors reproduce the independently printed ORCA R to the
    # precision allowed by its five-decimal vector components.
    vector_dot = (
        row.electric_x * row.magnetic_x
        + row.electric_y * row.magnetic_y
        + row.electric_z * row.magnetic_z
    )
    assert math.isclose(vector_dot, -420.54720, rel_tol=5e-4)


def test_cli_writes_orca_cpl_csv(tmp_path):
    output_path = tmp_path / "orca.out"
    output_path.write_text(ORCA_OUTPUT, encoding="utf-8")
    csv_path = tmp_path / "orca-cpl.csv"

    main(
        [
            "--log",
            str(output_path),
            "--program",
            "auto",
            "--cpl_only",
            "--cpl_csv_path",
            str(csv_path),
        ]
    )

    with csv_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 2
    assert rows[0]["state"] == "1"
    assert math.isclose(float(rows[0]["g value"]), build_cpl_rows(ORCA_OUTPUT)[0].g_value)
