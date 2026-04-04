# Code Review: DBS Paired Pulse R Analysis Scripts

## Overview

This document summarizes the statistical audit and fixes applied to the paired-pulse DBS evoked potential R analysis pipeline. Three scripts were reviewed: the main multi-subject analysis, the single-subject anesthesia analysis (3d413), and the single-subject conditioning length analysis (a23ed).

All scripts now use `_PairedPulseData_avg_5.csv` (peak-to-peak, 5-trial average) and save output to `DBS_EP_PairedPulse/R_output/`.

---

## Main Script: `dose_response_R_script_trim_conditions.R`

### Data Structure

```
12 subjects (9 with 1 channel, 3 with 2 channels)
  -> 16 subject:channel combinations
     -> 60 blocks total (2-6 per subject:channel)
        -> 2 stim levels per block (after min_stim_level filter)
           -> ~28 trials per cell (after 5-trial averaging and trimming)
```

6,464 trial-level rows aggregate to 120 cell means.

### Primary Model: Reduced Block RE (Approach 3)

```r
fit.lmmPP.block.reduced = lmerTest::lmer(PPvec ~ mapStimLevel + chanInCond +
    overallBlockType*pre_post + (1|subjectNum:chanVec:blockVec), data=dataList)
```

- Trial-level data with block random effect (60 groups)
- Drops `disease` and `baLabel` (underpowered: 12 subjects for disease, 16 channels for 2-df baLabel)
- Subject/channel-level variation absorbed by block RE
- Satterthwaite df ~50 for condition effects (not 6,464)
- No singular fit

### Comparison Models

| Approach | Data | Random Effects | Obs | Notes |
|----------|------|---------------|-----|-------|
| 1: Cell-mean (full) | `dataListAgg` | `(1\|subjectNum/chanVec)` | 120 | Singular fit (channel variance ~0) |
| 2: Block RE (full) | `dataList` | `(1\|subjectNum:chanVec:blockVec)` | 6,464 | All covariates |
| 3: Block RE (reduced) | `dataList` | `(1\|subjectNum:chanVec:blockVec)` | 6,464 | **Primary** |

Side-by-side comparison saved to `R_output/tab_model_comparison.html`.

### emmeans and Effect Sizes

- `emm_options(lmer.df = "asymptotic")` set for the block RE model to avoid slow Satterthwaite on 6,464 obs
- All plot CIs use `t_crit = qt(0.975, edf_conservative)` where `edf_conservative = ngrps - length(fixef)` (~52 df)
- `emm_options` reset to `"satterthwaite"` before smaller auxiliary models (`fit.lmmdiff`, `fit.effectSize`)
- Effect sizes computed via `eff_size()` with conservative df

### Plots

1. Interaction plot: overallBlockType x pre_post (emmip)
2. EMM by condition with CIs
3. Pre-post contrast within each condition (95% CI)
4. Standardized effect sizes (Cohen's d) for pre-post contrasts (log scale)
5. Fixed-effect comparison across all three model approaches

---

## Single-Subject: `baseline_variability_3d413.R`

### Design

Single subject (3d413), 2 channels (4, 6), 5 blocks per channel (3 asleep, 2 awake), ~18 trials per cell after 5-trial averaging.

### Model

```r
fit.lmm1 = lmerTest::lmer(PPvec ~ mapStimLevel + blockType + chanVec + (1|blockVec), data=dataList)
```

- Trial-level data with block RE (10 groups -- sufficient)
- `chanVec` as fixed effect (only 2 channels, too few for RE)
- `blockType` (awake/asleep) is the effect of interest

### Effect Size Analysis

Cohen's d (awake vs asleep) computed per (channel, stim_level) on trial-level data. Quantifies effect magnitude relative to trial-to-trial measurement noise.

Note: PPvec is log-transformed (`log_data=TRUE`), so d is on log scale.

---

## Single-Subject: `length_conditioning_a23ed.R`

### Design

Single subject (a23ed), 1 channel (5), 4 blocks with **separate baselines**:
- Comparison 1: block 2 (baseline) vs block 3 (post 5-min A/B 200 ms)
- Comparison 2: block 5 (baseline) vs block 6 (post 15-min A/B 200 ms)

The two baselines are different measurement blocks separated in time.

### Models

**Model 1 -- Full interaction:**
```r
fit.lmmPP.full = lmerTest::lmer(PPvec ~ mapStimLevel + overallBlockType*pre_post + (1|blockVec), data=dataList)
```
Singular fit expected: 4 fixed-effect combinations saturate 4 blocks. Kept for emmeans comparison.

**Model 2 -- Reduced (`pre_post` only with block RE):**
```r
fit.lmmPP.prepost = lmerTest::lmer(PPvec ~ mapStimLevel + pre_post + (1|blockVec), data=dataList)
```
Tests whether conditioning changes EP overall (pooling across durations). Non-singular but only 2 df for `pre_post`.

Side-by-side comparison saved to `R_output/tab_model_a23ed_comparison.html`.

### Why `blockVec` Cannot Be a Fixed Effect

Each block maps to exactly one combination of `overallBlockType*pre_post`. Including both as fixed effects creates a rank-deficient design matrix. This is the fundamental constraint of the single-subject design with 4 blocks.

### Permutation Tests (Primary Statistical Inference)

Non-parametric permutation tests (10,000 shuffles) provide model-free inference by shuffling trial block labels within each stim level:

**Test 1 -- 5-min post vs baseline (per stim level):** Shuffles trials between blocks 2 and 3.

**Test 2 -- 15-min post vs baseline (per stim level):** Shuffles trials between blocks 5 and 6.

**Test 3 -- Interaction (15-min effect minus 5-min effect):** Independently shuffles within each comparison, then computes the difference of differences. Tests whether longer conditioning produces a larger effect.

Results:
- 5-min conditioning: marginal at stim 2 (p=0.06), not significant at stim 3 (p=0.28), significant at stim 4 (p=0.008)
- 15-min conditioning: highly significant at all stim levels (p < 0.001)
- 15-min minus 5-min: significant at all stim levels (p=0.003, p<0.001, p=0.016)

P-values are two-sided: `mean(abs(perm_diffs) >= abs(obs_diff))`.

### Effect Size Analysis

Cohen's d (post vs baseline) computed per (stim_level, conditioning_length) on trial-level data. PPvec is NOT log-transformed (`log_data=FALSE`), so d is on original uV scale.

### Plots

1. Per-subject scatter plots by block (existing)
2. Across-condition scatter and summary plots (existing)
3. EMM interaction plot (overallBlockType x pre_post)
4. **Cohen's d bar plot** with CI error bars and small/medium/large reference lines
5. **Permutation null distribution histogram: 5-min post vs baseline**
6. **Permutation null distribution histogram: 15-min post vs baseline**
7. **Permutation null distribution histogram: 15-min effect minus 5-min effect (interaction)**
8. **Two-panel subplot**: each comparison with its own baseline plotted separately (shared y-axis)

---

## Issues Fixed

| Issue | Severity | Status |
|-------|----------|--------|
| Trial-level pseudoreplication (main) | High | **Fixed** -- block RE |
| Trial-level pseudoreplication (3d413) | High | **Fixed** -- block RE |
| Trial-level pseudoreplication (a23ed) | High | **Fixed** -- block RE + permutation tests |
| AR(1) degenerate/overparameterized | Medium | **Fixed** -- switched to lmer |
| `fit.lmmdiff` missing channel nesting | Medium | **Fixed** |
| Covariates dropped from glmmTMB | Medium | **Fixed** -- intentional reduction |
| Wrong model in avgMeas residual plots | Medium | **Fixed** |
| `sidVec`/`diseaseVec` redefinitions | Low-Med | **Fixed** -- single definition |
| `emm_options` global side effect | Medium | **Fixed** -- scoped per model |
| Hardcoded Windows path (3d413) | Low | **Fixed** -- removed |
| Typo `blocktype` (3d413) | Low | **Fixed** |
| Duplicate contrast/eff_size calls (a23ed) | Low | **Fixed** |
| Latent bug: rmsVec on pk_pk file | Low | **Fixed** -- all scripts use PPvec |

## Known Limitations

- **a23ed singular fit**: The full interaction model has 4 blocks for 4 condition combinations. Block RE is inherently unestimable. Permutation tests and effect size plots are the primary analysis.
- **Cell-mean model singular fit**: `(1|chanVec:subjectNum)` estimated near-zero with only 3 multi-channel subjects. Not the primary model.
- **`edf_conservative` is ad hoc**: `ngrps - length(fixef)` is a reasonable between-within approximation but not a standard formula. For publication, Satterthwaite df on key contrasts would be more defensible.
- **Cohen's d at trial level**: Assumes independent trials within blocks. Temporal autocorrelation between trials may inflate d. Standard caveat for single-subject neuroscience effect sizes.
- **emmip CIs use asymptotic df** while other plots use `t_crit` with conservative df. Minor visual inconsistency.
- **Stim level ordering differs by data file**: The `_rms_pk_pk` data showed stim level 3 > stim level 4 for a23ed (possible saturation). The current `_avg_5` peak-to-peak data shows a monotonic dose-response (stim 4 > stim 3 > stim 2), consistent with expectations.

## Data Files

All scripts use `{sid}_PairedPulseData_avg_5.csv` (peak-to-peak amplitude, 5-trial sequential average). Column `PPvec` is in volts, converted to microvolts (*1e6) in the scripts.

### Trial Averaging (MATLAB)

The `_avg_5` files are generated by the MATLAB pipeline via `avg_every_p_elems.m` (in `matlab_ecog_code` or `betaOscillationTriggerStimPaper` repos). This function averages every `p` consecutive trials along the trial dimension:

1. Trials are grouped into consecutive windows of size `p` (e.g., trials 1-5, 6-10, ...)
2. Each group is averaged to produce one "averaged trial"
3. **Remainder handling**: if the total trial count is not evenly divisible by `p`, the leftover trials are averaged together as a final shorter group (e.g., with 31 trials and p=5: 6 groups of 5 + 1 group of 1)

This means the last averaged trial may have higher variance than others when the remainder is small (1-2 trials). With `p=5` and typical block sizes of 30 trials, this produces 6 averaged trials per block per stim level.

### CSV File Variants

| Suffix | Metric | Averaging | Columns |
|--------|--------|-----------|---------|
| `_new_rms_pk_pk.csv` | RMS | none (individual trials) | PPvec, rmsVec |
| `_new_rms_pk_pk_avg_3.csv` | RMS | 3-trial sequential | PPvec, rmsVec |
| `_new_pk_pk_avg_5.csv` | peak-to-peak | 5-trial sequential | PPvec |
| `_avg_5.csv` | peak-to-peak | 5-trial sequential | PPvec |

The `_avg_5` and `_new_pk_pk_avg_5` files have identical column structure (PPvec only, no rmsVec). All current scripts use `_avg_5.csv`.

## Output

All figures and HTML tables save to `DBS_EP_PairedPulse/R_output/` when `savePlot = 1`. Key outputs:

- `tab_model_comparison.html` -- main script: 3-model side-by-side table
- `tab_model_a23ed_comparison.html` -- a23ed: full vs reduced model table
- `prepost_effsize_AUC.png` -- main script: effect sizes by condition
- `model_comparison_AUC.png` -- main script: fixed-effect forest plot
- `subj_a23ed_perm_*.png` -- permutation null distribution histograms (5-min, 15-min, interaction)
- `subj_a23ed_two_panel_prepost.png` -- two-panel pre/post with separate baselines
- `subj_a23ed_effect_size_conditioning_length.png` -- Cohen's d bar plot
- `subj_3d413_effect_size_awake_asleep.png` -- Cohen's d bar plot (log scale)
