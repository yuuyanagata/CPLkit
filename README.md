# CPLkit

`CPLkit.py` is a Python script that generates transition electric dipole moment density cubes, transition magnetic dipole moment density cubes, and CPL summary data from Gaussian TD-DFT output files. The script reads Gaussian excited-state information, reconstructs TEDM and TMDM densities on a cube grid, and exports CPL-related quantities to a CSV file. The implementation uses NumPy for numerical operations and may generate molecular orbital cube files through Gaussian `cubegen`, or it may consume pre-generated MO cubes directly.

## Features

- Generation of TEDM density cubes from Gaussian TD-DFT output.
- Generation of TMDM density cubes from Gaussian TD-DFT output.
- Automatic export of CPL-related quantities for all parsed excited states to CSV.
- Automatic export of the radiative rate constant for each transition to the CSV file as `Radiative Rate Constant (ns^-1)`.
- Support for either checkpoint based MO cube generation through `cubegen` or reuse of pre-generated `mo<MO>.cube` files.
- Preservation of Gaussian cube metadata for compatibility with GaussView.
- Optional export of Cartesian component cubes in addition to projected total cubes.

## Requirements

### Python

- Python 3.10 or later is recommended.

### Python packages

Install the Python dependency with:

```bash
pip install -r requirements.txt
```

### External software

One of the following workflows is required.

1. Gaussian `cubegen` is available and a checkpoint or formatted checkpoint file is provided through `--chk`.
2. Pre-generated molecular orbital cube files are available in a directory specified by `--mocubes_dir`.

## Installation

Clone the repository and install the Python dependency.

```bash
git clone https://github.com/yuuyanagata/CPLkit
cd CPLkit
pip install -r requirements.txt
```

## Input

The script expects a Gaussian TD-DFT output file supplied through `--log`.

For cube generation, one of the following is additionally required.

- `--chk` for checkpoint based MO cube generation via `cubegen`
- `--mocubes_dir` for a directory that contains `mo<MO>.cube` files

## Basic usage

### 1. Export only the CPL CSV

```bash
python CPLkit.py --log sample.out --cpl_only
```

### 2. Generate TEDM and TMDM cubes for an excited state using `cubegen`

```bash
python CPLkit.py --log sample.out --state 1 --chk sample.fchk
```

### 3. Generate cubes from pre-generated MO cubes

```bash
python CPLkit.py --log sample.out --state 1 --mocubes_dir ./sample_S1_mocubes
```

### 4. Also write Cartesian component cubes

```bash
python CPLkit.py --log sample.out --state 1 --chk sample.fchk --keep_components
```

## Radiative rate constant in the CSV

The CPL CSV includes a column named `Radiative Rate Constant (ns^-1)` for each parsed excited-state transition. The value is estimated from the oscillator strength and transition wavelength according to

```math
k_{\mathrm r}(\mathrm{ns}^{-1}) = 6.67 \times 10^{4} \frac{f}{\lambda^{2}(\mathrm{nm})}
```

where `f` is the oscillator strength reported by Gaussian TD-DFT and `了` is the transition wavelength in nm written to the same CSV row. This expression corresponds to the Einstein spontaneous emission rate in a practical wavelength based form.

## Main command-line options

- `--log`  
  Gaussian TD-DFT output file.

- `--state`  
  Excited-state index for cube generation. This option is not required when `--cpl_only` is used.

- `--chk`  
  Checkpoint or formatted checkpoint file for `cubegen`.

- `--mocubes_dir`  
  Directory containing `mo<MO>.cube` files.

- `--cubegen`  
  Path to the `cubegen` executable. The default value is `cubegen`.

- `--cubegen_npts`  
  First numerical argument passed to `cubegen`.

- `--cubegen_grid`  
  Remaining grid arguments passed to `cubegen`. The default value is `-3 h`.

- `--outdir`  
  Output directory.

- `--outprefix`  
  Output prefix. The default value is `S<state>`.

- `--overwrite_mo_cubes`  
  Regenerate MO cubes even if cached cube files are already present.

- `--keep_components`  
  Write `x`, `y`, and `z` component cubes in addition to projected total cubes.

- `--dtype`  
  Internal floating-point precision. Available values are `float32` and `float64`.

- `--coords`  
  Coordinate handling mode. Available values are `auto`, `aligned`, and `general`.

- `--cpl_only`  
  Export only the CPL CSV and skip cube generation.

- `--no_cpl_csv`  
  Disable CPL CSV export.

- `--cpl_csv_path`  
  Explicit output path for the CPL CSV file.

## Output

Depending on the selected options, the script writes files such as the following.

- `<logstem>-CPL.csv` containing CPL summary quantities and the `Radiative Rate Constant (ns^-1)` column for each transition
- `<logstem>_S1_TEDM_total.cube`
- `<logstem>_S1_TMDM_total.cube`
- `<logstem>_S1_TEDM_x.cube`
- `<logstem>_S1_TEDM_y.cube`
- `<logstem>_S1_TEDM_z.cube`
- `<logstem>_S1_TMDM_x.cube`
- `<logstem>_S1_TMDM_y.cube`
- `<logstem>_S1_TMDM_z.cube`

When `cubegen` is used, the script also creates a molecular orbital cube cache directory of the form:

- `<logstem>_S1_mocubes/`

## Repository structure

```text
CPLkit.py
LICENSE
README.md
README_ja.md
requirements.txt
```

## Citation and attribution

If this script is reused or adapted, citation of the original repository and indication of any modifications are appreciated as a matter of scholarly practice. The MIT License does not require academic citation, although preservation of the copyright notice and license text is required in copies or substantial portions of the software.

## License

This repository is distributed under the MIT License.

See the [LICENSE](LICENSE) file for details.
