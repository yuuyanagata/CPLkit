# Validation mode and contribution-map mode

CPLKit separates two density conventions because they answer different questions.

## Validation mode

Validation mode uses coefficient-linear transition-density amplitudes. It is designed to test whether real-space integration of EDTM and MDTM density components reproduces the Gaussian transition-moment tables. This mode should be used when the output is intended as a numerical validation of the reconstruction procedure.

Default factors are:

```text
EDTM weight = -2*c
MDTM scale  = -1*c
```

The electric prefactor includes the electron charge sign in the EDTM density convention. The MDTM expression uses the antisymmetrized orbital-gradient form implemented in `cplkit.densities`.

## Contribution-map mode

Contribution-map mode uses signed square configuration weights. It is intended for visualization of signed spatial contributions from dominant one-electron transitions. The magnitude follows the closed-shell configuration contribution, and the sign follows the printed CI coefficient.

Default weights are:

```text
EDTM weight = 2*c*abs(c)
MDTM scale  = c*abs(c)
```

Contribution-map mode is not a coefficient-linear transition-density amplitude and should not be interpreted as a strict numerical reconstruction of Gaussian EDTM or MDTM vector components.
