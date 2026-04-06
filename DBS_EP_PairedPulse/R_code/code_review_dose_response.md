# Code Review: DBS Paired Pulse R Analysis Scripts

## Overview

This document summarizes the statistical audit and fixes applied to the paired-pulse DBS evoked potential R analysis pipeline. Three scripts were reviewed:

1. **Primary analysis** (multi-subject): `dose_response_R_script_trim_conditions.R` — 12 subjects, tests whether paired-pulse conditioning protocols change EP amplitude
2. **Single-subject anesthesia analysis**: `baseline_variability_3d413.R` — subject 3d413, awake vs asleep EP variability
3. **Single-subject conditioning length analysis**: `length_conditioning_a23ed.R` — subject a23ed, 5-min vs 15-min conditioning duration

All scripts use `_PairedPulseData_avg_5.csv` (peak-to-peak, 5-trial average) and save output to `DBS_EP_PairedPulse/R_output/`.

---

## 1. Primary Analysis: `dose_response_R_script_trim_conditions.R`

### Data Structure

```
12 subjects (9 with 1 channel, 3 with 2 channels)
  -> 16 subject:channel combinations
     -> 60 blocks total (2-6 per subject:channel)
        -> 2 stim levels per block (after min_stim_level filter)
           -> ~12 trials per cell (after 5-trial averaging, no trimming)
```

1,462 trial-level rows (after filtering to 4 conditions: A/A 200, A/A 25, A/B 200, A/B 25) aggregate to 124 cell medians. 37 unique subject:block recording sessions; 13 unique subject:channel combinations; 9 subjects. (More cells pass the n>=10 filter with `trim_data = FALSE` than with trimming.)

A/B 100 is excluded (only 1 subject). A/A 25 was previously excluded without documented rationale but is now included -- it is the natural same-polarity control for A/B 25 (matched ISI, 3 subjects).

| Condition | Subjects | Channels | Blocks | Cell medians |
|-----------|----------|----------|--------|-------------|
| A/A 200 | 2 | 2 | 4 | 8 |
| A/A 25 | 3 | 4 | 6 | 16 |
| A/B 200 | 8 | 12 | 16 | 44 |
| A/B 25 | 5 | 6 | 12 | 32 |

### Robustness Strategy: Median Without Trimming

All scripts use `trim_data = FALSE` and aggregate trial-level data to **cell medians** rather than cell means. The median is inherently robust to outliers without requiring trimming, and combining both (trim + median) would be redundant -- removing extreme trials before taking a statistic that already ignores extremes wastes data without additional benefit. With ~12 trials per cell, every trial matters for precision.

This approach is applied consistently across all three scripts (main, 3d413, a23ed) including the permutation tests in a23ed.

### Log Transformation

EP amplitude (PPvec) is log-transformed before cell-median aggregation (`log_data = TRUE`). Log transformation is required for two reasons:

1. **Variance stabilization for the crossed RE.** On the raw (uV) scale, the crossed RE model is **singular** -- the `subjectNum:blockVec` variance component collapses to exactly zero because EP amplitude exhibits a mean-variance relationship (larger EPs = more variance). This causes the block RE to be unestimable and inflates Satterthwaite df from ~35 to ~103, making all block-level tests anti-conservative. On the log scale, variance is stabilized and all three RE components (subject SD=0.43, channel SD=0.61, block SD=0.18) are estimable.

2. **Standard practice for EP amplitude data.** EP amplitudes are bounded at zero with occasional large responses, producing right-skewed distributions. Log is the standard transformation.

Residual diagnostics on 124 cell medians (current primary model, log scale): Shapiro-Wilk p=0.006, skewness=-0.21, excess kurtosis=1.19. The mild departure from normality is expected at n=124 and lmer is robust to this level.

### Primary Model: Cell-Median with Crossed RE

```r
fit.lmmPP = lmerTest::lmer(PPvec ~ mapStimLevel + chanInCond +
    overallBlockType*pre_post + (1|subjectNum/chanVec) + (1|subjectNum:blockVec),
    data=dataListAgg)
```

**Cell-median aggregation:** Trials within a cell (subject x channel x block x stim level) are not independent -- they are sequential recordings from the same neural state. Aggregating to cell medians (124 obs from 1,462 trials) eliminates this pseudoreplication entirely.

**Crossed RE structure:** Multi-channel subjects (6, 7, 11, 12) had both channels recorded simultaneously during the same blocks. Blocks are therefore **crossed** with channels, not nested within them. The RE structure is:

| RE term | Groups | What it captures |
|---------|--------|-----------------|
| `(1\|subjectNum)` | 9 | Between-subject differences |
| `(1\|subjectNum:chanVec)` | 13 | Channel-level variance (where chanInCond lives) |
| `(1\|subjectNum:blockVec)` | 37 | Block-level variance (shared recording session across channels) |

Variance components: subject SD=0.43, channel SD=0.61, block SD=0.18, residual SD=0.24. Not singular.

`disease` and `baLabel` dropped (underpowered: 9 subjects for disease, only 1 MD; 13 channels for 3-level baLabel).

**Why this model over alternatives:** The flat `(1|subjectNum:chanVec:blockVec)` RE on trial-level data gives anti-conservative df for channel-level effects (df=37 for chanInCond instead of the correct ~7) and doesn't account for simultaneously recorded channels sharing block-level state. The cell-median with channel-only RE `(1|subjectNum/chanVec)` omits the block RE, inflating significance for block-level effects (df=65 instead of ~35). The crossed cell-median model gives correct df at every level.

### Comparison Models

| Approach | Data | Random Effects | Obs | Notes |
|----------|------|---------------|-----|-------|
| **Cell-median crossed RE** | `dataListAgg` | `(1\|subj/chan) + (1\|subj:block)` | 124 | **Primary** |
| Cell-median chan only | `dataListAgg` | `(1\|subj/chan)` | 124 | Anti-conservative for block-level effects |
| Trial-level flat RE | `dataList` | `(1\|subj:chan:block)` | 1,462 | Anti-conservative for channel-level effects |

Side-by-side comparison saved to `R_output/tab_model_comparison.html`.

### ANOVA

Type III ANOVA with Satterthwaite df (`lmerTest::anova(fit, type=3)`). Type III is used because the model includes an interaction term (`overallBlockType*pre_post`); Type II assumes no interaction and would be inappropriate.

Key results from Type III ANOVA (primary model):

| Effect | F | df | p |
|--------|---|-----|---|
| mapStimLevel | 51.2 | 1, 78 | 4.0e-10 |
| chanInCond | 8.3 | 1, 7 | 0.024 |
| overallBlockType | 2.44 | 3, 39 | 0.079 |
| pre_post | 0.57 | 1, 33 | 0.46 |
| overallBlockType:pre_post | 2.59 | 3, 36 | 0.068 |

No significant conditioning effects. The interaction is marginal (p=0.068) but not robust across analytic choices (see Sensitivity Analysis below).

### Degrees of Freedom by Hierarchy Level

Satterthwaite df reflect where each fixed effect lives in the hierarchy:

| Level | Effect | df | Source |
|-------|--------|-----|--------|
| Residual (within-block) | mapStimLevel | ~78 | 124 obs - 37 block groups - params. 2 stim levels per block. Lowest level of hierarchy, no further nesting needed. |
| Block (within-subject) | overallBlockType, pre_post, interaction | ~33-39 | 37 unique recording sessions minus params. Conditioning protocol and pre/post status vary between blocks. |
| Channel (within-subject) | chanInCond | ~7 | 13 subject:channel combinations, only 3 subjects with both chanInCond levels. |

The crossed RE structure ensures Satterthwaite assigns appropriate df at each level. A flat `(1|subj:chan:block)` would give df~37 for all effects, which is anti-conservative for channel-level effects.

### Posthoc Tests

**Pre vs post within each condition (Kenward-Roger df):**

| Condition | Estimate (post - pre) | Cohen's d | df | p |
|-----------|----------------------|-----------|-----|---|
| A/A 200 | 0.07 | 0.09 | 30 | 0.68 |
| A/A 25 | 0.22 | 0.28 | 34 | 0.24 |
| A/B 200 | 0.17 | 0.21 | 27 | 0.12 |
| A/B 25 | -0.23 | -0.29 | 31 | 0.068 |

No significant pre-to-post changes in any condition. All effect size CIs cross zero. Directional trends: A/B 200 trending up (small d=0.21), A/B 25 trending down (small d=-0.29), A/A conditions near zero. The opposite-polarity (A/B) conditions show larger directional effects than same-polarity (A/A) controls, consistent with the hypothesis but not statistically reliable with this sample size.

### emmeans and Confidence Intervals

- `emm_options(lmer.df = "satterthwaite")` -- fast on 124 obs, no need for asymptotic approximation
- All plots (emmip, EMM by condition, pre-post contrast, effect sizes) use **emmeans-native 95% CIs** with Satterthwaite df per cell/contrast. This avoids a single approximate t_crit and instead uses the correct df for each specific estimate.
- The emmip interaction plot uses `position_dodge` to horizontally offset overlapping conditions for visual clarity.

### Effect Sizes (Cohen's d)

Effect sizes for pre-post contrasts are computed via `emmeans::eff_size()`.

**Denominator choice: total SD.** The denominator is the total SD reconstructed from all model variance components:

```r
sigma_total <- sqrt(sum(as.numeric(VarCorr(fit.lmmPP))) + sigma(fit.lmmPP)^2)
# = sqrt(subj_var + chan_var + block_var + resid_var)
```

Note: `as.numeric(VarCorr(...))` returns **variances** (not SDs) in lme4 -- verified empirically.

The total SD was chosen over the residual-only SD following the recommendation from the MRC Cognition and Brain Sciences Unit FAQ on effect sizes (imaging.mrc-cbu.cam.ac.uk/statswiki/FAQ/tdunpaired): classical Cohen's d using pooled/total SD is the only variant unaffected by experimental design, making it suitable for meta-analytic synthesis and the standard d benchmarks (small=0.2, medium=0.5, large=0.8).

d values are on the **log-transformed** EP magnitude scale since `log_data = TRUE`.

### Sensitivity Analysis

The conditioning interaction p-value is sensitive to analytic choices:

| Configuration | Interaction p |
|---------------|--------------|
| 3 conditions, trimmed, cell-mean, flat RE | 0.45 |
| 3 conditions, trimmed, cell-mean, crossed RE | 0.073 |
| 4 conditions, trimmed, cell-mean, crossed RE | 0.13 |
| 4 conditions, untrimmed, cell-median, crossed RE | **0.068** |

The interaction moves between 0.05 and 0.45 depending on the number of conditions, aggregation method, and RE structure. **This is not a robust finding.** The directional pattern (A/B conditions show larger pre-post effects than A/A controls) is consistent but does not reach conventional significance under any defensible model specification.

Robust findings (consistent across all configurations):
- **Stimulation level drives EP amplitude** (p < 1e-06 in every model)
- **Channels in the conditioning pair have larger EPs** (p ~ 0.02 in every model with correct df)
- **No individual conditioning protocol produces a significant pre-to-post change**

### Auxiliary Models

Two secondary models are also fit for supplementary analysis:

1. **`fit.lmmdiff`**: `absDiff ~ mapStimLevel + chanInCond + blockType + (1|subjectNum/chanVec)` -- models the absolute difference from baseline.

2. **`fit.effectSize`**: `effectSize ~ meanPP + mapStimLevel + chanInCond + blockType + (1|subjectNum)` -- regresses trial-level Cohen's d on covariates. Produces `boundary (singular) fit` warning.

QQ plots for all models are saved as `pairedPulse_qq_*_allSubjs_AUC.png`.

### Plots

All EMM-related plots are saved as both PNG (300 dpi) and EPS in square dimensions (7x7 in).

1. Interaction plot: overallBlockType x pre_post with 95% CIs (dodged)
2. EMM by condition with CIs
3. Pre-post contrast within each condition (95% CI)
4. Standardized effect sizes (Cohen's d, total SD) for pre-post contrasts
5. Fixed-effect comparison across all three model approaches

---

## 2. Single-Subject Anesthesia Analysis: `baseline_variability_3d413.R`

### Design

Single subject (3d413), 2 channels (4, 6), 5 blocks per channel (3 asleep, 2 awake), ~18 trials per cell after 5-trial averaging.

### Model

```r
fit.lmm1 = lmerTest::lmer(PPvec ~ mapStimLevel + blockType + chanVec + (1|blockVec), data=dataList)
```

- Trial-level data with block RE (10 groups -- sufficient for RE estimation)
- `chanVec` as fixed effect (only 2 channels, too few for RE)
- `blockType` (awake/asleep) is the effect of interest
- `trim_data = FALSE`, `log_data = FALSE`

Log transformation was tested but produces catastrophic residual non-normality for this subject (skewness=-5.7, kurtosis=65) due to near-zero EP values creating a long left tail in log space. Raw-scale residuals are much better behaved (skewness=-0.5, kurtosis=2.9). Log is needed for the multi-subject analysis (variance stabilization across subjects) but not for single-subject data with a more homogeneous variance structure.

### Effect Size Analysis

Cohen's d (awake vs asleep) computed per (channel, stim_level) on trial-level data. Raw uV scale.

---

## 3. Single-Subject Conditioning Length Analysis: `length_conditioning_a23ed.R`

### Design

Single subject (a23ed), 1 channel (5), 4 blocks with **separate baselines**:
- Comparison 1: block 2 (baseline) vs block 3 (post 5-min A/B 200 ms)
- Comparison 2: block 5 (baseline) vs block 6 (post 15-min A/B 200 ms)

The two baselines are different measurement blocks separated in time. `trim_data = FALSE`, `log_data = FALSE` (raw uV for interpretability).

### Models

The mixed model approach is fundamentally limited here: with only 4 blocks, there are too few groups for a random effect (need >= 5-6). Block cannot be a fixed effect because it is perfectly confounded with `overallBlockType*pre_post`. The mixed models in the script (full interaction and reduced pre_post-only) are included for completeness but have convergence warnings and should be interpreted with caution.

Side-by-side comparison saved to `R_output/tab_model_a23ed_comparison.html`.

### Permutation Tests (Primary Statistical Inference)

Non-parametric permutation tests (10,000 shuffles) provide the primary inference by shuffling trial block labels within each stim level. The test statistic is the **difference of medians**.

**Test 1 -- 5-min post vs baseline (per stim level):** Shuffles trials between blocks 2 and 3.

**Test 2 -- 15-min post vs baseline (per stim level):** Shuffles trials between blocks 5 and 6.

**Test 3 -- Interaction (15-min effect minus 5-min effect):** Independently shuffles within each comparison, then computes the difference of differences.

Results (median-based, untrimmed):
- 5-min conditioning: marginal at stim 2 (p=0.055), not significant at stim 3 (p=0.75), significant at stim 4 (p=0.045)
- 15-min conditioning: highly significant at all stim levels (p<0.001, p<0.001, p=0.001)
- 15-min minus 5-min interaction: marginal at stim 2 (p=0.053), not significant at stim 3/4 (p=0.14, p=0.18)

P-values are two-sided: `mean(abs(perm_diffs) >= abs(obs_diff))`.

**Important caveats:**

1. **Trial exchangeability.** The permutation test shuffles individual trials between blocks, assuming they are exchangeable under the null. This is the same independence assumption as a t-test or lm -- if trials within a block are temporally correlated, the permutation null distribution is too narrow and p-values may be anti-conservative. The 5-trial averaging upstream partially mitigates temporal autocorrelation but does not eliminate it.

2. **Permutation space.** With ~10 trials per block per stim level, the number of unique permutations is C(~20, ~10) ~ 184,756. The minimum achievable p-value is ~1/184,756 ~ 5e-06. The reported p<0.001 values are well above this floor.

3. **Per-stim-level vs pooled.** Per-stim-level tests (reported above) are the conceptually correct approach -- they never mix trials from different stimulation intensities. The pooled interaction test in the script mixes stim levels, which conflates stim-level composition with the conditioning effect. The per-stim-level interaction tests are underpowered (only ~10 trials per cell for the median).

### Effect Size Analysis

Cohen's d (post vs baseline) computed per (stim_level, conditioning_length) on trial-level data. PPvec is NOT log-transformed (`log_data=FALSE`), so d is on original uV scale. Note this makes effect sizes not directly comparable to the main multi-subject analysis (which uses log scale).

### Plots

1. Per-subject scatter plots by block
2. Across-condition scatter and summary plots
3. EMM interaction plot (overallBlockType x pre_post)
4. Cohen's d bar plot with CI error bars and small/medium/large reference lines
5. Permutation null distribution histogram: 5-min post vs baseline
6. Permutation null distribution histogram: 15-min post vs baseline
7. Permutation null distribution histogram: 15-min effect minus 5-min effect (interaction)
8. Two-panel subplot: each comparison with its own baseline (shared y-axis)

---

## Issues Fixed

| Issue | Severity | Status |
|-------|----------|--------|
| Trial-level pseudoreplication (main) | High | **Fixed** -- cell-median aggregation |
| Flat RE giving wrong df for channel effects | High | **Fixed** -- crossed RE structure |
| Simultaneously recorded channels not modeled | High | **Fixed** -- `(1\|subjectNum:blockVec)` |
| A/A 25 condition excluded without rationale | Medium | **Fixed** -- included (3 subjects) |
| Trial-level pseudoreplication (3d413) | High | **Fixed** -- block RE |
| Trial-level pseudoreplication (a23ed) | High | **Fixed** -- permutation tests |
| Trimming + median redundancy | Medium | **Fixed** -- `trim_data = FALSE` with median |
| Effect size using residual-only SD | Medium | **Fixed** -- total SD for comparability |
| CIs using single approximate t_crit | Low-Med | **Fixed** -- emmeans-native CIs per contrast |
| AR(1) degenerate/overparameterized | Medium | **Fixed** -- switched to lmer |
| `emm_options` global side effect | Medium | **Fixed** -- scoped per model |
| `as_data_frame()` deprecation | Low | Noted -- should migrate to `as_tibble()` |

## Known Limitations

- **Marginal interaction is not robust.** The conditioning x pre_post interaction (p=0.068) is sensitive to the number of conditions included, aggregation method, and RE structure. It should be reported as exploratory, not confirmatory.
- **a23ed: 4 blocks is insufficient for mixed models.** Block RE requires >= 5-6 groups. All mixed models for a23ed have convergence warnings. Permutation tests are the primary analysis, with the caveat that they assume trial exchangeability.
- **a23ed permutation tests assume trial independence.** Shuffling individual trials between blocks is the same exchangeability assumption as a t-test. Within-block temporal correlation (partially mitigated by 5-trial averaging) could make p-values anti-conservative.
- **Cohen's d on different scales.** Main script uses log(uV), a23ed uses raw uV. Effect sizes are not directly comparable between the two analyses.
- **Small sample size.** 9 subjects, 37 recording sessions, 124 cell medians. Power for block-level effects is limited by the ~37 block groups.
- **Disease and Brodmann area underpowered.** Only 1 MD patient (of 9) after filtering. baLabel has 3 levels across 13 channels. Neither can be reliably tested.
- **Unequal 5-trial averaging.** The last averaged trial in a block may be the average of 1-4 trials (when trial count is not divisible by 5), giving it higher variance. Cell-median aggregation partially mitigates this.
- **mapStimLevel as ordered factor.** Polynomial contrasts assume equally spaced levels, but actual stimulation currents (mA) are not equally spaced across subjects. The contrasts are meaningful for rank order only, not physical units.

## Data Files

All scripts use `{sid}_PairedPulseData_avg_5.csv` (peak-to-peak amplitude, 5-trial sequential average). Column `PPvec` is in volts, converted to microvolts (*1e6) in the scripts.

### Trial Averaging (MATLAB)

The `_avg_5` files are generated by the MATLAB pipeline via `avg_every_p_elems.m` (in `matlab_ecog_code` or `betaOscillationTriggerStimPaper` repos). This function averages every `p` consecutive trials along the trial dimension:

1. Trials are grouped into consecutive windows of size `p` (e.g., trials 1-5, 6-10, ...)
2. Each group is averaged to produce one "averaged trial"
3. **Remainder handling**: if the total trial count is not evenly divisible by `p`, the leftover trials are averaged together as a final shorter group (e.g., with 31 trials and p=5: 6 groups of 5 + 1 group of 1)

### CSV File Variants

| Suffix | Metric | Averaging | Columns |
|--------|--------|-----------|---------|
| `_new_rms_pk_pk.csv` | RMS | none (individual trials) | PPvec, rmsVec |
| `_new_rms_pk_pk_avg_3.csv` | RMS | 3-trial sequential | PPvec, rmsVec |
| `_new_pk_pk_avg_5.csv` | peak-to-peak | 5-trial sequential | PPvec |
| `_avg_5.csv` | peak-to-peak | 5-trial sequential | PPvec |

All current scripts use `_avg_5.csv`.

## Output

All figures and HTML tables save to `DBS_EP_PairedPulse/R_output/` when `savePlot = 1`. Key outputs:

- `tab_model_comparison.html` -- main script: 3-model side-by-side table
- `tab_model_a23ed_comparison.html` -- a23ed: full vs reduced model table
- `emmip_AUC.png/eps` -- main script: EMM interaction plot with CIs
- `prepost_effsize_AUC.png/eps` -- main script: effect sizes by condition
- `model_comparison_AUC.png/eps` -- main script: fixed-effect forest plot
- `subj_a23ed_perm_*.png/eps` -- permutation null distribution histograms
- `subj_a23ed_two_panel_prepost.png/eps` -- two-panel pre/post with separate baselines
- `subj_a23ed_effect_size_conditioning_length.png/eps` -- Cohen's d bar plot
- `subj_3d413_effect_size_awake_asleep.png/eps` -- Cohen's d bar plot (log scale)
