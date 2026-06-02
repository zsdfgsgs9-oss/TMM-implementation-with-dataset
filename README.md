# Transfer Matrix Method for Electron Beam Focusing

MATLAB implementation of the transfer matrix method (TMM) for simulating and optimizing magnetic-field-focused electron beam optical systems.

## Overview

This code simulates electron trajectories from a photocathode to a silicon wafer target in an electron beam lithography system. The transfer matrix method models the complete electron optical system as a series of 3×3 homogeneous coordinate transfer matrices, achieving ~2000× speedup compared to COMSOL FEM particle tracing.

## File Structure

### Main Scripts
- `main_process.m` — Main simulation entry point (spot mode)
- `optimization_U.m` — Accelerating voltage optimization
- `analyze_stripes_11_16.m` — Stripe pattern post-processing (TMM vs FEM comparison)
- `verify_space_charge.m` — Space charge model validation
- `plot_performance_comparison.m` — TMM vs FEM performance benchmark
- `convert_comsol_data.m` — COMSOL CSV to MATLAB format conversion

### Transfer Matrix Functions
Each function computes a 3×3×N transfer matrix for a specific physical effect:
- `calc_M_spread.m` — Energy spread & angular divergence
- `calc_M_EB0.m` — Zero-order E/B field error (defocus)
- `calc_M_Br1.m` — Radial B-field gradient (rotation distortion)
- `calc_M_E1.m` — Higher-order E-field error (defocus)
- `calc_M_B1.m` — B-field variation (image scaling)
- `calc_M_EB_angle.m` — E×B drift (stretch distortion)
- `calc_M_ee.m` — Space charge (electron-electron interaction)

### Utilities
- `init_particles.m` — Particle initial state generator
- `calc_common_physics.m` — Common physics pre-computation
- `export_figs_final.m` — Figure export (600 DPI PNG/TIFF/EPS)
- `ee_calculation_exact.m` — Exact space charge integral calculation

### Data Files
- `*粒子初始状态.csv` — COMSOL-exported initial particle states
- `*焦面图像.csv` — COMSOL-exported target plane images
- `*.fig` — MATLAB figure files for paper

## Usage

```matlab
% Main simulation
main_process

% Stripe pattern analysis (files 11-16)
analyze_stripes_11_16

% Export publication-quality figures
export_figs_final
```

## Requirements
- MATLAB R2020a or later (uses `pagemtimes` for batched 3D matrix multiplication)
- Signal Processing Toolbox (for `findpeaks`)
- Curve Fitting Toolbox (optional, for `fit`)

## Reference
See `电子束聚焦公式推导（新）.docx` for the complete physics derivation.
