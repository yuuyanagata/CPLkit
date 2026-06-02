@echo off
setlocal

REM Replace these paths with local files.
set "LOG=path\to\calculation.out"
set "MOCUBES=path\to\mocubes"
set "OUTDIR=results\contribution_map"

python ..\CPLkit.py ^
  --log "%LOG%" ^
  --state 1 ^
  --mocubes_dir "%MOCUBES%" ^
  --outdir "%OUTDIR%" ^
  --outprefix S1 ^
  --density_mode contribution-map ^
  --cube_units bohr ^
  --accumulation vectorized ^
  --vectorized_stack_mode memmap ^
  --keep_components ^
  --skip_validation_report

pause
