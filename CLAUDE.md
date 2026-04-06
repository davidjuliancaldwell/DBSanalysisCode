# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

MATLAB and R code to analyze evoked potentials (EPs) from intraoperative Deep Brain Stimulation (DBS) experiments and paired pulse stimulation experiments in humans. The pipeline extracts EP waveforms from raw recordings in MATLAB, then performs statistical analysis (linear mixed models) in R.

Author: David J Caldwell. BSD-3 License.

## Data Pipeline

```
Raw TDT Recording → MATLAB Load/Convert → Stimulus Extraction →
Optional ICA Artifact Removal → Peak-to-Peak Detection →
Optional Trial Averaging (avg_every_p_elems) →
CSV Export → R Statistical Analysis → Figures
```

## Key Entry Points

- **`DBS_EP_PairedPulse/master_script_analyze_EP.m`** — Main paired pulse EP analysis. Set `sidVecIterate` to choose subjects, configure flags (`savePlot`, `saveData`, `screenBadChans`, `tryArtifact`), then run.
- **`DBS_EP_PairedPulse/prepare_EP_blocks.m`** — Per-subject parameter config (large switch on `sid`). Defines `stimChans`, `tBegin/tEnd`, `badTrials`, `blockLabel` for each subject.
- **`DBS_EP_PairedPulse/R_code/dose_response_R_script_trim_conditions.R`** — Main R statistical analysis across all subjects (N=9). Cell-median aggregation with crossed RE `(1|subjectNum/chanVec) + (1|subjectNum:blockVec)` on log-transformed EP amplitude. 4 conditioning protocols (A/A 200, A/A 25, A/B 200, A/B 25). See `R_code/code_review_dose_response.md` for full audit and `R_code/statistical_results_summary.md` for results with all numbers.
- **`DBS_EP_PairedPulse/R_code/baseline_variability_3d413.R`** — Single-subject anesthesia variability (awake vs asleep). Block RE model + trial-level Cohen's d effect sizes. `trim_data = FALSE`, `log_data = FALSE` (log produces catastrophic non-normality for this subject).
- **`DBS_EP_PairedPulse/R_code/length_conditioning_a23ed.R`** — Single-subject conditioning length (5 vs 15 min). Permutation tests with median statistic (primary) + effect sizes. `trim_data = FALSE`, `log_data = FALSE`.

## Architecture

### Constants / Configuration (3 tiers)
- `Z_ConstantsDBS.m` — Internal EP subjects, global `OUTPUT_DIR` and `META_DIR`
- `DBS_EP_PairedPulse/Z_ConstantsDBS_PairedPulse.m` — Paired pulse subjects, dirs
- `externalEPs/Z_ConstantsDBS_externalEPs.m` — External electrode subjects

All constants files depend on environment variables `dbs_subject_dir` and `OUTPUT_DIR` (via `myGetenv`).

### Analysis Functions (`analysisFunctions/`)
- `peak_to_peak.m` — Core peak-to-peak amplitude algorithm
- `extract_PP*.m` — Three wrappers: averaged, single-trial, index-from-average
- `extract_rms_single_trial.m` — RMS amplitude analysis
- `plot_EPs.m` — EP waveform visualization
- `savitzkyGolay.m` — Savitzky-Golay filter implementation

### Trial Averaging
`avg_every_p_elems.m` (external, in `matlab_ecog_code` or `betaOscillationTriggerStimPaper` repos) averages consecutive trials in groups of `p`. Remainder trials (when count not divisible by `p`) are averaged as a final shorter group. The `_avg_5` CSV files use `p=5`.

### Subject Identification
Subjects use hashed 5-character IDs (e.g., `c963f`, `46c2a`, `3d413`). ~20 paired pulse patients, ~10 external EP patients. Adding a new subject requires entries in the relevant `Z_Constants` file and a new case in `prepare_EP_blocks.m`.

### Data Directory Layout
Subjects are stored under `dbs_subject_dir` with structure:
```
SubjectID/MATLAB_Converted/EP_Measurement/[block_data].mat
```

Output CSVs go to `DBS_EP_PairedPulse/R_data/` (~128 files). R analysis output (figures, HTML tables) goes to `DBS_EP_PairedPulse/R_output/`. Per-subject MNI electrode coordinates live in `DBS_EP_PairedPulse/R_config_files/`.

### Other Pipelines
- **`internal_DBS_EP/`** — Opposite-polarity artifact subtraction for internal electrodes
- **`externalEPs/`** — Simpler analysis for external surface electrodes
- **`DBS_EP_PairedPulse/intraoperative/`** — Real-time EP visualization during surgery
- **`ica_train_dbs.m` / `ica_artifact_remove_train_dbs.m`** — FastICA-based artifact removal (root level)
- **`DBS_EP_PairedPulse/vizualization/`** — MNI coordinate and response plots using FreeSurfer surfaces

## R Statistical Analysis Notes

- **Primary model**: Cell-median aggregation (no trimming) with crossed RE on log-transformed EP amplitude. See `R_code/code_review_dose_response.md` for full audit and `R_code/statistical_results_summary.md` for results.
- **Robustness strategy**: `trim_data = FALSE`, median aggregation across all scripts. `log_data = TRUE` (main only, required for crossed RE convergence) / `FALSE` (3d413 and a23ed, single-subject raw scale).
- **Key documentation files**:
  - `R_code/code_review_dose_response.md` — Statistical audit, model specifications, known limitations
  - `R_code/statistical_results_summary.md` — All numbers for reviewers: ANOVA tables, effect sizes, CIs, sensitivity analysis

## R Dependencies

plyr, here, nlme, ggplot2, drc, minpack.lm, lmtest, glmm, lme4, multcomp, lmerTest, sjPlot, emmeans, dplyr, effectsize, caret, DescTools, report, easystats

R scripts use `here()` for path resolution from the repo root.

## MATLAB Dependencies

Signal Processing Toolbox, Statistics Toolbox. FastICA (external) required only if `tryArtifact = 1`. TDT data loading via `promptForTDTrecording.m`.

## Git / File Conventions

- `.gitignore` excludes: `*.png *.mat *.asv *.zip *.fig *.svg *.eps *.pdf *.Rhistory *.RData`
- FreeSurfer imaging data lives in `DBS_EP_PairedPulse/imaging/` (13 patients)
- Third-party helpers: `cbrewer_helper/` (color palettes), `jbfill/` (confidence interval shading)
