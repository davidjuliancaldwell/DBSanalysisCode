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
- **`DBS_EP_PairedPulse/R_code/dose_response_R_script_trim_conditions.R`** — Main R statistical analysis across all subjects (N=13). Cell-median aggregation with crossed RE + random linear dose slope per channel `(1+stim_linpoly|subjectNum:chanVec) + (1|subjectNum) + (1|subjectNum:blockVec)` on log-transformed EP amplitude. 4 conditioning protocols (A/A 200, A/A 25, A/B 200, A/B 25). 262 cell medians, 17 channels. Amplitude filter: 10-1000 uV. All stim levels included (`min_stim_level=1`); all subjects mapped to ordered levels 1-4. Also writes per-subject, per-condition modulation CSVs (`modulation_all_per_subj_cond.csv`, `modulation_depression_levels34.csv`, `modulation_augmentation_levels34.csv`) to `R_output/`; values are on raw uV scale (cell medians exponentiated from log scale since `log_data = TRUE`). See `R_code/code_review_dose_response.md` for full audit and `R_code/statistical_results_summary.md` for results with all numbers.
- **`DBS_EP_PairedPulse/R_code/baseline_variability_3d413.R`** — Single-subject anesthesia variability (awake vs asleep). Trial-level lmer with block RE + uncorrelated random linear dose slope `(1|blockVec) + (0+stim_linpoly|blockVec)`. Permutation tests (primary). `trim_data = FALSE`, `log_data = FALSE` (log produces catastrophic non-normality for this subject). Amplitude filter 10-1000 uV, all stim levels.
- **`DBS_EP_PairedPulse/R_code/length_conditioning_a23ed.R`** — Single-subject conditioning length (5 vs 15 min). Permutation tests with median statistic (primary) + effect sizes. `trim_data = FALSE`, `log_data = FALSE`. Per-stim-level p values adjusted with Benjamini-Hochberg FDR within each of three families (5-min, 15-min, interaction; 4 stim levels each); `perm_p_BH` column appears alongside `perm_p` in the console output and the docx tables. Stratified pooled tests remain uncorrected (single omnibus test per family).
- **`DBS_EP_PairedPulse/R_code/permutation_tests_all_subjects.R`** — Applies the same per-stim-level trial-level median-difference permutation test (10,000 iterations, stratified by stim level, two-sided) to every (subject, channel, pre/post pair) in `R_config_files/subj_*.R`, with BH-FDR correction within each (subject, channel, condition) family across stim levels. Uses the same sidVec as the main analysis (excludes 3d413). Writes `R_output/permutation_results_all_subjects.csv` (all rows) and `R_output/permutation_results_all_subjects_stim4.csv` (max-stim subset).

## Architecture

### Constants / Configuration (3 tiers)
- `Z_ConstantsDBS.m` — Internal EP subjects, global `OUTPUT_DIR` and `META_DIR`
- `DBS_EP_PairedPulse/Z_ConstantsDBS_PairedPulse.m` — Paired pulse subjects, dirs
- `externalEPs/Z_ConstantsDBS_externalEPs.m` — External electrode subjects

All constants files depend on environment variables `dbs_subject_dir` and `OUTPUT_DIR` (via `myGetenv`).

### Environment Setup
- **`DBS_EP_PairedPulse/setupEnvironment.m`** — Sets `dbs_subject_dir` and `OUTPUT_DIR` environment variables and adds required paths (`MATLAB_ECoG_code`, this repo, and `helpers/` last so local copies shadow same-named externals). Run once per MATLAB session before any analysis scripts, or call at the top of a script.
- **`helpers/`** — Local copies of external functions patched for cross-platform use. Added to the MATLAB path last so they take precedence (`addpath` prepends).
  - `SaveFig.m` — Patched to accept Unix absolute paths (and any Windows drive letter). Original in `MATLAB_ECoG_code/Visualization/SaveFig.m` only recognized `C:/` or `D:/` and silently wrote to `c:/Tim/research/script/generated_figs/<full_unix_path>/` on macOS, creating a literal `c:` directory in the cwd. External copy also patched.
  - `TouchDir.m` — Local copy of the dependency used by `SaveFig`.

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
- **`DBS_EP_PairedPulse/visualize_EP_by_stimlevel.m`** — Mean EP waveforms by stim level for specified subjects (pooled across blocks, per-block, and with CI shading)
- **`DBS_EP_PairedPulse/visualize_PP_extraction.m`** — Visualize peak-to-peak extraction pipeline (avg 5 → SG smooth → peak_to_peak) with peak/trough markers. Configurable subject list; includes switch with all 14 subjects' block/channel configs matching R analysis
- **`DBS_EP_PairedPulse/visualize_EP_pre_post_conditioning.m`** — Pre vs post conditioning EP waveforms at configurable stim levels. USER CONFIG block at top: `stimInterest` (default `[3 4]`), `cases` (default = top 3 depression cases from `modulation_depression_levels34.csv`: 9f852 chan 4 with A/B 25 ms, 68574 chan 7 with A/A 200 ms, 41a73 chan 5 with A/A 200 ms; chan 7 of 9f852 excluded from the pipeline due to atypical morphology), `subjectBlocks`. Per case produces one figure with `2 × length(stimInterest)` tiled panels (row 1 = mean only, row 2 = mean + 95% CI; pre = grey, post = orange). Caches `prepare_EP_blocks` output per subject.
- **`DBS_EP_PairedPulse/visualize_conditions_comparison.m`** — Reproduces the 2019 `a23ed_chan_5_mean_measure.fig` style for any list of `(sid, channel, blocks, legendText)` cases. Loops over `cases` (caching `prepare_EP_blocks` output per subject) and, at each stim level in `stimInterest`, emits a mean-only figure and a 95% CI figure. Mean figure carries heavy dashed `vline`s at the baseline-block peak and trough (from `pkLocsBlock{1}` / `trLocsBlock{1}` at that stim level; CI figure has no vlines). Uses cbrewer `Dark2` colors. Saves `.png` (600 dpi), `.eps` (vector), and `.fig` with filenames `{sid}_chan_{ch}_blocks_{b1}_{b2}_..._{mean|confInt}_measure_stim{s}.*`. Default style: square 7×7 in figure, `pbaspect([1 1 1])`, YLim `[-300 300]`, XTicks `[0 20 40 60]`, YTicks `[-200 -100 0 100 200]`, FontSize 22 with title multiplier 0.8. Per-subject YLim override via `yLimPerSid` / `yTicksPerSid` Maps (e.g., `68574 → [-140 140]` with yticks `[-100 0 100]` because its EP amplitudes are much smaller than other subjects). Default `cases`: a23ed `[2 3 6 8]` (all-conditions overlay) and `[5 6]` (matched pre-post pair used in the R statistical test for 15-min A/B 200 ms), 41a73 `[5 6]`, 68574 `[6 7]`, 9f852 `[2 3]`.

## R Statistical Analysis Notes

- **Primary model**: Cell-median aggregation (n>=3 per cell, no trimming) with crossed RE + channel dose slope on log-transformed EP amplitude. 13 subjects (9 PD, 4 MD; includes a23ed 15-min sessions), **230 cell medians, 16 channels** (was 262/17 before 9f852 chan 7 was dropped on 2026-04-20). `mapStimLevel` is an ordered factor (1-4) for all subjects; previous special mapping for 41a73/68574 (`c(1,3,4,0)`) removed — all 4 stim levels now included with standard ordering. Data regenerated from committed MATLAB code (latest 2026-04-20). See `R_code/code_review_dose_response.md` for full audit and `R_code/statistical_results_summary.md` for results.
- **Contrast coding**: `mapStimLevel` uses `contr.poly` (ordered, centered). Unordered factors in interaction terms (`overallBlockType`, `pre_post`, `blockType`, `halfBlock`) use `contr.sum` (sum-to-zero), set explicitly per-variable before model fitting — ensures Type III main effects test averaged over the other factor's levels. Additive-only factors and RE grouping factors use R default `contr.treatment`. No global contrast options set.
- **Robustness strategy**: `trim_data = FALSE`, median aggregation across all scripts. `log_data = TRUE` (main only, required for crossed RE convergence) / `FALSE` (3d413 and a23ed, single-subject raw scale). Amplitude filter 10-1000 uV applied in R (MATLAB pipeline does not filter by magnitude).
- **MATLAB data pipeline**: per-trial per-channel baseline subtraction ([-50 −5] ms pre-stim window, applied to `epochsEP` in `prepare_EP_blocks.m:875`) → waveform averaging (5 trials) → windowing (tBegin-tEnd) → Savitzky-Golay smoothing (order 3, frame 91) → peak-to-peak extraction. PP amplitude is `|max_pos| + |max_neg|` within `[tBegin, tEnd]`, so it is **invariant** to any constant per-trial DC shift — baseline subtraction affects plotted waveforms but leaves `PPvec` unchanged (narrowing `[tBegin, tEnd]` itself, however, does change `PPvec` because it restricts which extrema `findpeaks` can see). The [-50 -5] baseline window matches the 2019 reference figures and the beta-osc paper convention (`B_ExtractNeuralData_PP_reref.m:309`). R uses `PPvec` column. Data regenerated via `DBS_EP_PairedPulse/regenerate_avg5_data.m` (latest: 2026-04-20).
- **Per-subject extraction windows (2026-04-20):** a23ed `tBegin=3, tEnd=30` ms (was 5/50 — wider window was catching a late slow rebound as the "peak"); 9f852 `tBegin=3, tEnd=30` ms (was 7/70 — same issue, plus the old tBegin=7 missed the early positive deflection). Other subjects unchanged.
- **9f852 channel 7 excluded from the pipeline (2026-04-20):** morphology was not consistent with a typical EP (flat, late extrema); removed from `R_config_files/subj_9f852.R` (`chanIntVec = c(4)`), both viz scripts (`visualize_EP_pre_post_conditioning.m`, `visualize_conditions_comparison.m`), and `vizualization/mni_DBS_electrodeLoc_plots.m`. Drops the primary analysis from 17 subject-channels (262 cells) to 16 subject-channels (230 cells).
- **Key documentation files**:
  - `R_code/code_review_dose_response.md` — Statistical audit, model specifications, known limitations
  - `R_code/statistical_results_summary.md` — All numbers for reviewers: ANOVA tables, effect sizes, CIs, sensitivity analysis
  - `R_output/statistical_tables.docx` / `statistical_tables_3d413.docx` / `statistical_tables_a23ed.docx` — Manuscript-ready tables (extracted from models)
  - `R_output/descriptive_stats_baseline_diff.docx` — Raw uV and percent differences from baseline
  - `R_code/descriptive_stats_docx.R` — Script to regenerate descriptive stats
  - `R_output/main_analysis_workspace.RData` — Saved R workspace with all model objects for quick loading
  - `R_output/modulation_all_per_subj_cond.csv` — Per-subject, per-condition, per-channel, per-stim-level cell-median pre/post with absolute and percent differences and direction label (`augmentation` / `depression` / `no_change`). Generated by `dose_response_R_script_trim_conditions.R`.
  - `R_output/modulation_depression_levels34.csv` / `modulation_augmentation_levels34.csv` — Subsets of the above filtered to stim levels 3-4 with negative / positive modulation, sorted by magnitude.

## R Dependencies

plyr, here, nlme, ggplot2, drc, minpack.lm, lmtest, glmm, lme4, multcomp, lmerTest, sjPlot, emmeans, dplyr, effectsize, caret, DescTools, report, easystats

R scripts use `here()` for path resolution from the repo root.

## MATLAB Dependencies

Signal Processing Toolbox, Statistics Toolbox. FastICA (external) required only if `tryArtifact = 1`. TDT data loading via `promptForTDTrecording.m`.

## Git / File Conventions

- `.gitignore` excludes: `*.png *.mat *.asv *.zip *.fig *.svg *.eps *.pdf *.Rhistory *.RData`
- FreeSurfer imaging data lives in `DBS_EP_PairedPulse/imaging/` (13 patients)
- Third-party helpers: `cbrewer_helper/` (color palettes), `jbfill/` (confidence interval shading), `helpers/` (local patched copies of `SaveFig` + `TouchDir` for cross-platform path handling)
