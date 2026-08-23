# CPLKit

CPLKit is an independently written Python package for circularly polarized luminescence oriented post-processing of Gaussian and ORCA excited-state calculations. It exports CPL summary tables from Gaussian or ORCA TD-DFT/TDA output. For Gaussian calculations it can additionally reuse or generate molecular-orbital cube files and reconstruct electric dipole transition moment (EDTM) and magnetic dipole transition moment (MDTM) density grids.

CPLKit was implemented de novo in Python from published transition-density equations and Gaussian/ORCA output conventions. It is not a MATLAB port, wrapper, or line-by-line translation.

## Features

- Gaussian TD-DFT excited-state parser for CI coefficients and transition-moment tables.
- CPL summary CSV export from Gaussian or ORCA TD-DFT/TDA output.
- Automatic Gaussian/ORCA output detection and ORCA length-gauge spectrum parsing.
- Native Gaussian `cubegen` orchestration for molecular-orbital cube generation.
- Explicit Bohr handling for Gaussian cube files.
- EDTM and MDTM density cube generation.
- Separate `validation` and `contribution-map` density modes.
- Component-level validation against Gaussian EDTM and MDTM transition-moment tables.
- Low-memory streaming accumulation and faster NumPy vectorized accumulation.

## Requirements

CPLKit requires Python 3.9 or later and NumPy. Gaussian and `cubegen` are not distributed with CPLKit. They are needed only when molecular-orbital cube files must be generated from checkpoint files. If the required MO cube files already exist, CPLKit can reuse them through `--mocubes_dir`.

## Installation

Install from a local checkout.

```bash
python -m pip install .
```

After installation, the command-line entry point is available as:

```bash
cplkit --help
```

The repository also includes backward-compatible script entry points:

```bash
python CPLkit.py --help
```

## ORCA CPL summary

For an ORCA TD-DFT or TDA output containing the length-gauge absorption and CD
spectrum tables, export all matching excited states with:

```bash
cplkit \
  --log path/to/orca_calculation.out \
  --program orca \
  --cpl_only \
  --outdir results/orca_cpl
```

The default `--program auto` detects ORCA automatically. ORCA support is for
CPL CSV export; Gaussian-only MO cube reconstruction is not attempted. See
[docs/ORCA_CPL.md](docs/ORCA_CPL.md) for the required output tables, unit
convention, and the `g_CPL` definition.

## Density modes

CPLKit intentionally separates two density modes.

`validation` mode uses coefficient-linear transition-density amplitudes. It is intended for quantitative comparison between real-space EDTM and MDTM integrals and the corresponding Gaussian transition-moment tables.

`contribution-map` mode uses signed square configuration weights. The EDTM weight is `2*c*abs(c)`, and the MDTM scale is `c*abs(c)`. This convention follows the closed-shell configuration contribution `2*c**2` while preserving the CI-coefficient sign. It is intended for signed spatial contribution maps and is not expected to reproduce Gaussian transition-moment components by direct integration.

## Basic usage with existing MO cube files

The following command generates signed contribution maps from an existing Gaussian output file and a directory that already contains `mo<MO>.cube` files.

```bash
cplkit \
  --log path/to/calculation.out \
  --state 1 \
  --mocubes_dir path/to/mocubes \
  --outdir results/contribution_map \
  --outprefix S1 \
  --density_mode contribution-map \
  --cube_units bohr \
  --accumulation vectorized \
  --vectorized_stack_mode memmap \
  --keep_components \
  --skip_validation_report
```

The corresponding validation-mode calculation is:

```bash
cplkit \
  --log path/to/calculation.out \
  --state 1 \
  --mocubes_dir path/to/mocubes \
  --outdir results/validation \
  --outprefix S1 \
  --density_mode validation \
  --phase_align edtm \
  --cube_units bohr \
  --accumulation vectorized \
  --vectorized_stack_mode memmap \
  --keep_components
```

The principal visualization outputs are:

```text
*_EDTM_total.cube
*_MDTM_total.cube
```

When `--keep_components` is used, CPLKit also writes the x, y, and z component cubes for each transition moment.

## Generating MO cubes with cubegen

When checkpoint files are available, CPLKit can call Gaussian `cubegen`.

```bash
cplkit \
  --log path/to/calculation.out \
  --state 1 \
  --chk path/to/calculation.chk \
  --cubegen cubegen \
  --cubegen_grid "-3 h" \
  --outdir results/from_chk \
  --density_mode validation \
  --cube_units bohr
```

The Gaussian IOp(9/40) coefficient-printing threshold should be set sufficiently tightly for density reconstruction. The validation study associated with this code used IOp(9/40=4) for the reported Gaussian TD-DFT cases.

## Repository layout

```text
CPLkit.py                  Backward-compatible script entry point
cplkit/                    Python package
examples/                  Minimal command-line examples
docs/                      Implementation and mode notes
tests/                     Lightweight import and weighting tests
pyproject.toml             Packaging metadata
LICENSE                    MIT License
```

## License

CPLKit is distributed under the MIT License. Gaussian, Gaussian `cubegen`, and any third-party quantum-chemistry data files are not included in this repository and are not licensed by this repository.
