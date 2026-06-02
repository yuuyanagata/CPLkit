# CPLkit.py implementation notes

## Scope

CPLkit.py is a CPL-oriented Gaussian TD-DFT post-processing package. It parses Gaussian excited-state output, generates or reuses molecular-orbital cube files, reconstructs electric and magnetic transition-moment-density components, writes projected density cube files, and exports CPL summary CSV files.

## Code provenance

The present software was written de novo in Python from published transition-density equations and Gaussian output conventions. It is not a MATLAB port, wrapper, or line-by-line translation. No MATLAB source code is required, invoked, or distributed by this package.

## What is new in the present package

The novelty of this software is the integrated validated workflow rather than the known formula for gCPL. The package adds:

1. a modular Gaussian parser for transition-moment tables and excited-state configurations,
2. native Gaussian cubegen orchestration with concurrent molecular-orbital cube generation,
3. explicit Bohr unit handling for Gaussian cube files,
4. coefficient-printing diagnostics for IOp(9/40),
5. separate coefficient-linear validation and signed square contribution-map modes,
6. full TD-DFT response-density convention diagnostics based on X+Y and X-Y combinations,
7. component-level validation against Gaussian EDTM and MDTM tables,
8. CSV outputs that preserve vector-level, component-level, and best-convention validation results.

## Relation to general parsers and wavefunction-analysis programs

CPLkit.py is not intended to replace general tools such as cclib or Multiwfn. Those programs address broad parsing or wavefunction-analysis tasks. CPLkit.py addresses a narrower CPL-specific workflow in which Gaussian TD-DFT coefficient output, molecular-orbital cube files, response-density reconstruction, and transition-moment validation must be handled consistently.

## Density-mode rationale

Contribution-map mode uses signed square weights, `2*c*abs(c)`, because the magnitude follows the closed-shell configuration contribution `2*c**2` while the sign follows the printed CI coefficient. This convention is useful for displaying signed spatial contributions from dominant orbital transitions. It is not a coefficient-linear transition-density amplitude.

Validation mode uses coefficient-linear amplitudes because Gaussian EDTM and MDTM vector components are coherent transition matrix elements. Real-space integration is expected to reproduce the Gaussian transition-moment tables only in validation mode.
