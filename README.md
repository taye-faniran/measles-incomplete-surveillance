# Measles Transmission Inference Under Incomplete Surveillance

This repository contains the verified MATLAB reproducibility materials for:

> Faniran, T. S., Blumberg, S., and Cardenas, S. D. **Recorded Measles Cluster Sizes Inform Subcritical Transmission Only Jointly with Surveillance.** Under review at *BMC Medical Research Methodology*.

The repository includes the final analysis programs, numerical results, saved analysis objects, manuscript figures, grouped U.S. measles data, computational logs, and software-environment information used for the submitted study.

## Scientific question

Recorded outbreak clusters can be smaller than the underlying transmission clusters because surveillance misses cases. Consequently, recorded cluster sizes do not identify transmission intensity and individual-level heterogeneity independently of the observation process.

The analyses link:

1. a negative-binomial branching process with effective reproduction number `R_eff` and dispersion parameter `k`; and
2. a two-stage observation process with passive-detection probability `p_p` and active-detection probability `p_a` after surveillance is triggered within a cluster.

The code evaluates known-surveillance likelihood inference, fixed-data sensitivity to assumed surveillance, joint likelihood estimation with auxiliary surveillance information, MAP estimation with detection priors, and an application to grouped U.S. measles cluster summaries.

## What is included

| Folder | Contents |
|---|---|
| `MATLAB_Code` | Complete analysis and figure-generation programs |
| `Data_and_Numerical_Results` | Parameter estimates, profile surfaces, posterior samples, validation summaries, and model checks |
| `Saved_Analysis_Objects` | MATLAB checkpoints and verified saved results for computationally intensive analyses |
| `Tables` | Grouped U.S. measles data and final parameter/prior tables |
| `Figures` | Final manuscript figures in PDF and PNG format; MATLAB FIG files are included where available |
| `Documentation` | Analysis logs, MCMC settings, diagnostics, and the original MATLAB software environment |

`Reproducibility_Manifest_Original_BMC_Release.csv` records the contents of the submitted release. `SHA256_Original_BMC_Release.csv` preserves the checksums from that release. `PUBLIC_REPOSITORY_MANIFEST.csv` and `SHA256_Public_Repository.csv` record the contents and checksums of this portable public package.

## Software requirements

The submitted analyses were run in MATLAB R2020a on 64-bit Windows. The main required products are:

- MATLAB;
- Optimization Toolbox; and
- Statistics and Machine Learning Toolbox.

Parallel Computing Toolbox is optional for programs that can use parallel execution. See `Documentation/MATLAB_Software_Environment.txt` for the complete submitted environment.

## Quick start

Clone or download the repository, open MATLAB, and set the repository root as the current folder:

```matlab
addpath(fullfile(pwd,'MATLAB_Code'))
```

All newly generated material is written beneath an automatically created `outputs` folder. The public programs use repository-relative paths; no user-specific folder editing is required.

### Inspect the completed study without recomputation

The fastest audit path is:

1. inspect the final figures in `Figures`;
2. inspect parameter estimates and sensitivity results in `Data_and_Numerical_Results` and `Tables`;
3. inspect the analysis logs in `Documentation`; and
4. inspect the saved MATLAB objects in `Saved_Analysis_Objects`.

This path does not repeat the long simulations, profile calculations, or MCMC runs.

### Redraw Figure 1 from verified saved results

```matlab
Figure1_Known_Surveillance_Analysis
```

This program loads the packaged Figure 1 analysis object and CSV cross-check, then redraws the known-surveillance figure without resimulation or re-estimation.

### Recompute the fixed-data and JLE analyses

```matlab
Figure2_Fixed_Data_Analysis
Figure3_JLE_16Panel_Correct
```

Figure 2 generates one fixed recorded dataset and analyzes it repeatedly under different assumed detection settings. Figure 3 performs the joint likelihood analysis using auxiliary surveillance information.

After Figure 3 finishes, the Figure 4 workflow can be run in order:

```matlab
Repair_JLE_Global_Optimum_Then_Run_Figure4
Figure4_JLE_Detection_Profiles_Final
Figure4_JLE_Detection_Profiles_ManuscriptFinal
```

These programs audit the Figure 3 optimum and construct the coarse, refined, and final detection-profile results.

### Recompute the MAP analyses

```matlab
Run_MAP_Analysis
Redraw_MAP_Figures_5_7
Figure6_MAP_Detection_Profiles_ManuscriptFinal
```

The workflow is checkpointed. If a long computation stops, rerunning the same program reuses compatible completed checkpoints and continues unfinished work.

### Recompute the grouped U.S. analysis

```matlab
Run_US_Grouped_MAP_Analysis
```

The primary application uses completed years 2020–2023. The 2024 snapshot is treated as a prespecified sensitivity analysis because an outbreak was ongoing at the publication cutoff. The program uses grouped summaries and does not invent unavailable individual cluster sizes.

## Parameters

| Manuscript symbol | Common code name | Meaning |
|---|---|---|
| `R_eff` | `R`, `R0`, or `R_eff` | Effective reproduction number under the subcritical branching-process model |
| `k` | `k` | Negative-binomial dispersion; smaller values represent greater individual-level transmission heterogeneity |
| `p_p` | `p_p` or `p1` | Passive case-detection probability |
| `p_a` | `p_a` or `p2` | Active detection probability among cases missed passively after a cluster is triggered |

## Reproducibility design

- Simulation programs set explicit random-number seeds.
- Computationally intensive workflows save checkpoints.
- Profile likelihood or profile posterior calculations use direct numerical evaluation rather than cosmetic contour smoothing.
- Saved CSV results provide machine-readable numerical evidence behind the figures.
- The U.S. analysis includes simulation validation, prior sensitivity, posterior summaries, posterior-predictive checks, MCMC settings, and diagnostics.
- The original submitted release contained 89 checksum-listed files; all 89 were present and matched their recorded SHA-256 hashes before the public path-only portability changes were made.

## Computational cost

Several analyses are intentionally intensive. Full reproduction may require hours or longer depending on hardware, especially the multistart optimization, adaptive profile calculations, simulation validation, and MCMC workflows. Use the saved objects and numerical outputs when the objective is to audit the completed study rather than repeat every calculation.

## License

The code is released under the MIT License. See `LICENSE`.

## Contact

Taye Samuel Faniran, PhD  
University of California, San Francisco  
Taye.Faniran@ucsf.edu
