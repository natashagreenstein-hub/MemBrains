# Metabolism Matters: Scaling Noise by Regional Oxygen Consumption in Whole-Brain Models

**Paper:** Greenstein N, Firaux T, Aquilue D, Epp S, Castrillon Guzman G, Riedl V, Tuszynski J, Cavaglia M, Deco G. "Metabolism matters: scaling noise by regional oxygen consumption in whole-brain models improves fit to empirical data."

**Lead contact:** natasha.greenstein@gmail.com

---

## Overview

This repository contains the MATLAB code and processed data to reproduce all results in the paper. We show that scaling the per-region noise parameter β of a Hopf whole-brain model by the regional cerebral metabolic rate of oxygen (CMRO2) significantly improves model fit to empirical resting-state functional connectivity, and that this improvement depends on the specific anatomical arrangement of the metabolic map.

The model is fitted to group-averaged resting-state data from the Human Connectome Project (HCP, N=498) at the DBS80 parcellation (N=80 regions). Metabolic maps (CMRO2, CMRglc, OEF) are derived from quantitative fMRI in Epp et al. (2026, *Nature Neuroscience*).

---

## Repository Structure

```
code/
  NEMO_oxglc_AABB.m         Main analysis: fits β = β₀(B + A·mᵢ) for each
                             metabolic map (ox/glc ratio, CMRO2, CMRglc) via
                             PSO; runs uniform-shuffle and spatial-
                             autocorrelation-preserving (Moran) null models;
                             bootstraps error bars. Produces all three figures.

  hopf_int.m                Hopf oscillator integrator used during the RS-GEC
                             fitting loop (Lyapunov/analytic FC computation).

  hopf_int_pert_beta.m      Computes model FC analytically (Lyapunov equation)
                             given a per-region β vector. Called for every PSO
                             objective evaluation and null-model iteration.

  build_DBS80_geometry.m    Computes DBS80 parcel centroids in MNI coordinates
                             from the atlas NIfTI; output saved as
                             DBS80_geometry.mat and required by the spatial-
                             autocorrelation null in NEMO_oxglc_AABB.m.

  plot_NEMO_oxglc_results.m Loads results_NEMO_*_DBS80.mat files and produces
                             Figures 1–3.

data/
  SC_dbs80HARDIFULL.mat         DBS80 HARDI structural connectivity (80×80)
  empirical_HCP_rest_DBS80.mat  Pre-computed group FC, COV(τ), peak frequencies
  metabolic_maps_DBS80.mat      CMRO2, CMRglc, OEF per DBS80 region
  CMRO2_over_CMRglc_DBS80.mat  CMRO2/CMRglc ratio per region
  ox_gluc_ratio.mat             Normalised ox/glc ratio (unit mean, 1×80)
  DBS80_geometry.mat            Parcel centroids in MNI coordinates (80×3)
  results_NEMO_oxglc_DBS80.mat  PSO results: ox/glc ratio model (Fig. 1)
  results_NEMO_ox_DBS80.mat     PSO results: CMRO2-only model (Fig. 2)
  results_NEMO_glc_DBS80.mat    PSO results: CMRglc-only model (Fig. 2)
```

---

## Large Data Files (download separately)

The HCP resting-state timeseries is too large to include (~173 MB):

| File | Source |
|------|--------|
| `HCP_rest_DBS80_Functional_Timeseries.mat` | Human Connectome Project — http://www.humanconnectome.org/study/hcp-young-adult |

Download the HCP resting-state fMRI data (N=498 participants, 1200 volumes, TR=0.72 s, standard HCP minimal preprocessing pipeline) and place this file in the same directory as the code before running.

The metabolic map raw data (quantitative fMRI) is available at OpenNeuro: https://openneuro.org/datasets/ds004873. The parcellated versions are included in `data/` so re-parcellation is not required to reproduce the results.

---

## Requirements

- MATLAB R2026a or later
- No additional toolboxes required. The Moran spectral randomisation null is implemented directly in `NEMO_oxglc_AABB.m`.

---

## Reproducing the Results

Set `projDir` at the top of each script to point to this repository's root, then:

**Run the full analysis (all three metabolic maps → Figs. 1–3):**
```matlab
% Set metabolic_map before running — one of: 'ox_gluc', 'ox', 'glc'
metabolic_map = 'ox_gluc';   % Fig. 1
NEMO_oxglc_AABB
% Repeat for 'ox' and 'glc' (Fig. 2)
```

**Plot results:**
```matlab
plot_NEMO_oxglc_results
```

**Rebuild DBS80 geometry** (only needed if DBS80_geometry.mat is missing):
```matlab
build_DBS80_geometry
```

---

## Key Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| N | 80 | DBS80 parcels |
| β (baseline) | 0.01 | Homogeneous noise amplitude |
| a | −0.01 | Bifurcation parameter (sub-critical) |
| Bandpass | 0.008–0.08 Hz | 2nd-order Butterworth |
| TR | 0.72 s | HCP repetition time |
| τ | 3 TR | Lag for directed COV(τ) |
| PSO swarm size | 40 | Seeded near (B, A) = (0, 1), bounds [0, 5] |
| PSO iterations | 300 max | |
| Permutation nulls | 200 | Both uniform-shuffle and Moran spatial null |
| Bootstrap resamples | 200 | Subject resampling for error bars |

---

## Citation

If you use this code or data, please cite:

> Greenstein N, Firaux T, Aquilue D, Epp S, Castrillon Guzman G, Riedl V, Tuszynski J, Cavaglia M, Deco G. Metabolism matters: scaling noise by regional oxygen consumption in whole-brain models improves fit to empirical data. *[Journal]*, 2026. DOI: [TBD]

Please also cite the source datasets:

> Epp SM, Castrillón G, et al. BOLD signal changes can oppose oxygen metabolism across the human cortex. *Nature Neuroscience* 29, 1225–1236 (2026). https://doi.org/10.1038/s41593-025-02132-9

> Van Essen DC et al. The WU-Minn Human Connectome Project: an overview. *NeuroImage* 80, 62–79 (2013).

---

## License

Code released under the MIT License. HCP-derived data is subject to HCP Data Use Terms: https://www.humanconnectome.org/study/hcp-young-adult/data-use-terms
