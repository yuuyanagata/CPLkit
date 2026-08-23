import math

import cplkit
from cplkit.densities import excitation_weights


def test_version_is_public_semver():
    assert cplkit.__version__ == "0.5.0"


def test_contribution_map_weighting():
    edtm_weight, mdtm_scale = excitation_weights(0.5, "contribution-map")
    assert math.isclose(edtm_weight, 0.5)
    assert math.isclose(mdtm_scale, 0.25)


def test_validation_weighting_defaults():
    edtm_weight, mdtm_scale = excitation_weights(0.5, "validation")
    assert math.isclose(edtm_weight, -1.0)
    assert math.isclose(mdtm_scale, -0.5)
