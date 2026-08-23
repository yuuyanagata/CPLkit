# ORCA TD-DFT/TDA CPL summaries

CPLkit can export state-resolved CPL summary CSV files from ORCA TD-DFT and TDA
outputs. This path reads ORCA's final length-gauge absorption and circular
dichroism (CD) spectrum tables. It does not reconstruct real-space EDTM/MDTM
density cubes; that workflow remains specific to Gaussian `cubegen` data.

## Required ORCA output

The output must contain both of these tables:

```text
ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS
CD SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS
```

ORCA normally prints them for a TD-DFT/TDA spectrum calculation when CD output
is enabled. A representative input block is:

```text
%tddft
  TDA true
  NRoots 6
  DoCD true
end
```

Run CPLkit in summary-only mode:

```bash
cplkit --log calculation.out --program orca --cpl_only --outdir results
```

`--program auto` is the default, so `--program orca` can be omitted. When an
ORCA optimization prints several spectra, CPLkit uses the last absorption table
and the last CD table.

## Definition and units

For each excited state, CPLkit converts ORCA's electric transition dipole
components using

```text
mu_scaled = 254.17 * mu_ORCA(a.u.)
```

ORCA's printed length-gauge table convention obeys, within its output rounding,

```text
R(1e-40 cgs) = mu_scaled dot [2 * (-0.927401) * m_ORCA(a.u.)]
```

CPLkit therefore uses the separately printed rotatory strength `R` as the
signed scalar product and evaluates

```text
g_CPL = 4 * R / (|mu_scaled|^2 + |m_scaled|^2)
```

where `m_scaled = 2 * (-0.927401) * m_ORCA(a.u.)`. Using the printed `R`
preserves its sign and avoids additional rounding error from the vector
components. The CSV also reports the converted vectors, their norms,
`cos(theta)`, excitation energy, wavelength, oscillator strength, and the
radiative-rate estimate already used by the Gaussian summary path.

Only matching states present in both length-gauge tables are exported. Velocity
gauge tables are deliberately ignored so that gauges are never mixed.
