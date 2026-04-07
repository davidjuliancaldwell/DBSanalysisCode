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
12 subjects (8 with 1 channel, 4 with 2 channels)
  -> 16 subject:channel combinations
     -> 47 blocks total (after correcting 41a73 block labels)
        -> all stim levels included (min_stim_level = 1)
           -> n>=3 trials per (stim level, block) cell required
           -> amplitude filter: 10 < PPvec < 1000 uV
```

250 cell medians after filtering to 4 conditions (A/A 200, A/A 25, A/B 200, A/B 25), applying 10/1000 uV amplitude filter, and per-stim-level n>=3 cell filter. 51 unique subject:block recording sessions; 17 unique subject:channel combinations; 13 subjects (includes a23ed with 15-min conditioning sessions only; MNI coordinates unavailable, baLabel set to "Unknown"). Previously 9 subjects with n>=10 filter; relaxed to n>=3 with per-stim-level dropping (failing stim levels dropped individually rather than excluding entire subjects).

**Data regeneration (2026-04-07):** All `_PairedPulseData_avg_5.csv` files regenerated from committed MATLAB code (`regenerate_avg5_data.m`). Original files renamed to `_unknown_provenance`. Pipeline: waveform averaging (5 trials) -> windowing (tBegin-tEnd) -> Savitzky-Golay smoothing (order 3, frame 91) -> peak-to-peak extraction. No baseline normalization on PPvec. R uses the `PPvec` column (not `PPfromAvgVec`). Row counts differ from originals due to bad trial exclusions added in commit c55c8cc (April 2021).

**Amplitude filter:** Peak-to-peak values below 10 uV or above 1000 uV are discarded before analysis to exclude incorrectly extracted evoked potentials. The MATLAB extraction pipeline does not apply magnitude cutoffs; filtering is done in R after loading CSVs.

**Stimulation levels:** All available stimulation levels are included (`min_stim_level = 1`). The mixed model handles the resulting unbalanced design (different subjects contribute different numbers of stim levels).

A/B 100 is excluded (only 1 subject). A/A 25 was previously excluded without documented rationale but is now included -- it is the natural same-polarity control for A/B 25 (matched ISI, 3 subjects).

**41a73 block correction (2026-04-06):** Subject 41a73's R config file had block labels that disagreed with the current MATLAB `prepare_EP_blocks.m` definitions. The R config was created from an earlier version of the MATLAB block assignments and was never updated. Corrected to match MATLAB: only A/B 200 ms (block 3) and A/A 200 ms (block 6) conditioning protocols, with baselines at blocks 2 and 5 respectively. Previously, the R config incorrectly assigned A/B 25 and A/A 25 protocols to blocks that were baselines in the MATLAB definition.

### Robustness Strategy: Median Without Trimming

All scripts use `trim_data = FALSE` and aggregate trial-level data to **cell medians** rather than cell means. The median is inherently robust to outliers without requiring trimming, and combining both (trim + median) would be redundant -- removing extreme trials before taking a statistic that already ignores extremes wastes data without additional benefit. With ~12 trials per cell, every trial matters for precision.

This approach is applied consistently across all three scripts (main, 3d413, a23ed) including the permutation tests in a23ed.

### Log Transformation

EP amplitude (PPvec) is log-transformed before cell-median aggregation (`log_data = TRUE`). Log transformation is required for two reasons:

1. **Variance stabilization for the crossed RE.** On the raw (uV) scale, the crossed RE model is **singular** -- the `subjectNum:blockVec` variance component collapses to exactly zero because EP amplitude exhibits a mean-variance relationship (larger EPs = more variance). This causes the block RE to be unestimable and inflates Satterthwaite df from ~35 to ~103, making all block-level tests anti-conservative. On the log scale, variance is stabilized and all three RE components (subject SD=0.41, channel SD=0.63, block SD=0.21) are estimable.

2. **Standard practice for EP amplitude data.** EP amplitudes are bounded at zero with occasional large responses, producing right-skewed distributions. Log is the standard transformation.

### Primary Model: Cell-Median with Crossed RE + Channel Dose Slope

```r
fit.lmmPP = lmerTest::lmer(PPvec ~ mapStimLevel + chanInCond +
    overallBlockType*pre_post +
    (1 + stim_linpoly|subjectNum:chanVec) + (1|subjectNum) + (1|subjectNum:blockVec),
    data=dataListAgg)
```

**Cell-median aggregation:** Trials within a cell (subject x channel x block x stim level) are not independent -- they are sequential recordings from the same neural state. Aggregating to cell medians (232 obs, 12 subjects) eliminates this pseudoreplication entirely. Per-stim-level filter requires >=3 trials per (stim level, block) cell; stim levels that fail are dropped individually rather than excluding the entire subject.

**Crossed RE structure:** Multi-channel subjects (6, 7, 11, 12) had both channels recorded simultaneously during the same blocks. Blocks are therefore **crossed** with channels, not nested within them.

**Random slope on channel:** `stim_linpoly` is the linear polynomial contrast extracted from `contr.poly(n_stim_levels)[, 1]`, provided as a numeric column so it can be used as a random slope. Each electrode has its own dose-response steepness (slope SD=0.55), reflecting cortical position. Model comparison showed channel slope (AIC=158) beats block slope (AIC=210) and both-slopes (AIC=161, block slope collapses to SD~0.04 when channel slope is present). Dose-response steepness is a property of electrode location, not recording session.

| RE term | Groups | What it captures |
|---------|--------|-----------------|
| `(1\|subjectNum)` | 13 | Between-subject differences |
| `(1 + stim_linpoly\|subjectNum:chanVec)` | 17 | Channel-level intercept + dose-response slope |
| `(1\|subjectNum:blockVec)` | 51 | Block-level variance (shared recording session) |

Variance components: subject SD=0.33, channel intercept SD=0.60, channel slope SD=0.78 (corr=-0.11), block SD=0.22, residual SD=0.23. Not singular. Residuals pass Shapiro-Wilk normality (W=0.993, p=0.259, skewness=-0.22, excess kurtosis=0.24).

`disease` and `baLabel` dropped (underpowered: 9 subjects for disease, only 1 MD; 13 channels for 3-level baLabel).

**Why this model over alternatives:** The flat `(1|subjectNum:chanVec:blockVec)` RE on trial-level data gives anti-conservative df for channel-level effects (df=37 for chanInCond instead of the correct ~7) and doesn't account for simultaneously recorded channels sharing block-level state. The cell-median with channel-only RE `(1|subjectNum/chanVec)` omits the block RE, inflating significance for block-level effects (df=65 instead of ~35). The crossed cell-median model gives correct df at every level.

### Comparison Models

| Approach | Data | Random Effects | Obs | AIC | Notes |
|----------|------|---------------|-----|-----|-------|
| **Cell-median, channel slope** | `dataListAgg` | `(1+slope\|subj:chan) + (1\|subj) + (1\|subj:block)` | 202 | 158 | **Primary** |
| Cell-median, block slope | `dataListAgg` | `(1+slope\|subj:block) + (1\|subj/chan)` | 202 | 210 | Block slope collapses when channel slope present |
| Cell-median, intercept-only | `dataListAgg` | `(1\|subj/chan) + (1\|subj:block)` | 202 | 246 | No random slope |
| Trial-level flat RE | `dataList` | `(1\|subj:chan:block)` | 2,392 | 2384 | Anti-conservative for channel-level effects |

Side-by-side comparison saved to `R_output/tab_model_comparison.html`.

### ANOVA

Type III ANOVA with Satterthwaite df (`lmerTest::anova(fit, type=3)`). Type III is used because the model includes an interaction term (`overallBlockType*pre_post`); Type II assumes no interaction and would be inappropriate.

Key results from Type III ANOVA (primary model, channel dose slope):

| Effect | F | NumDF | DenDF | p |
|--------|---|-------|-------|---|
| mapStimLevel | 17.8 | 3 | 35.5 | 3.2e-07 |
| chanInCond | 8.25 | 1 | 9.9 | 0.017 |
| overallBlockType | 2.95 | 3 | 41.8 | 0.044 |
| pre_post | 1.27 | 1 | 35.0 | 0.268 |
| overallBlockType:pre_post | 1.42 | 3 | 38.5 | 0.251 |

Conditioning protocol main effect reaches significance (p=0.044). The A/B 200 vs A/B 25 Tukey pairwise comparison is significant (p=0.030). No significant interaction (p=0.251). Note: `mapStimLevel` has 3 NumDF because all stim levels are included (up to 4 levels per subject). The random channel slope absorbs between-electrode variability in dose-response steepness.

### Degrees of Freedom by Hierarchy Level

Satterthwaite df reflect where each fixed effect lives in the hierarchy:

| Level | Effect | df | Source |
|-------|--------|-----|--------|
| Channel (random slope) | mapStimLevel (.L) | ~12 | Linear dose-response slope varies across 13 channels |
| Residual (within-block) | mapStimLevel (.Q, .C) | ~147-150 | 202 obs - 41 block groups - params. Quadratic/cubic not absorbed by linear slope. |
| Block (within-subject) | overallBlockType, pre_post, interaction | ~27-33 | 41 unique recording sessions minus params. Conditioning protocol and pre/post status vary between blocks. |
| Channel (within-subject) | chanInCond | ~7 | 13 subject:channel combinations, only 3 subjects with both chanInCond levels. |

The crossed RE structure ensures Satterthwaite assigns appropriate df at each level. A flat `(1|subj:chan:block)` would give df~37 for all effects, which is anti-conservative for channel-level effects.

### Posthoc Tests

**Pre vs post within each condition (Kenward-Roger df):**

| Condition | Estimate (post - pre) | SE | 95% CI | df | Cohen's d | d 95% CI |
|-----------|----------------------|-----|---------|-----|-----------|----------|
| A/A 200 | -0.003 | 0.156 | [-0.319, 0.314] | 35.6 | -0.002 | [-0.301, 0.297] |
| A/A 25 | 0.203 | 0.186 | [-0.172, 0.577] | 41.8 | 0.19 | [-0.170, 0.545] |
| A/B 200 | 0.217 | 0.109 | [-0.004, 0.438] | 33.6 | 0.20 | [-0.018, 0.418] |
| A/B 25 | -0.078 | 0.121 | [-0.323, 0.167] | 37.7 | -0.07 | [-0.307, 0.163] |

No significant pre-to-post changes in any individual condition. All effect size CIs cross zero. A/B 200 approaches significance (p=0.054 from CI). Directional trends: A/B 200 trending up (small d=0.20), A/B 25 trending slightly negative (d=-0.07), A/A conditions near zero. Total SD used for d: 1.083.

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
| 4 conditions, untrimmed, cell-median, crossed RE (int only) | 0.068 |
| 4 conditions, untrimmed, cell-median, channel slope, n>=3, 13 subj, regenerated data (current) | **0.251** |

The interaction moves between 0.05 and 0.45 depending on the number of conditions, aggregation method, and RE structure. **This is not a robust finding.** The directional pattern (A/B conditions show larger pre-post effects than A/A controls) is consistent but does not reach conventional significance under any defensible model specification.

Robust findings (consistent across all configurations):
- **Stimulation level drives EP amplitude** (p < 1e-06 in every model)
- **Channels in the conditioning pair have larger EPs** (p ~ 0.02 in every model with correct df)
- **Conditioning protocol main effect** (p=0.044, A/B 200 vs A/B 25 Tukey p=0.030) — but confounded with between-subject differences since different subjects contribute to different protocols
- **No individual conditioning protocol produces a significant pre-to-post change** (A/B 200 approaches significance, d=0.20)
- **No significant conditioning x pre/post interaction** — the pre-to-post change does not significantly differ across protocols

### Secondary Model: halfBlock (Temporal Modulation)

Tests whether the conditioning effect differs between the first and second half of trials within each block (`fit.lmmPP.half`). Trials split at the midpoint within each (block, stim level) cell. Uses `dataListAggHalf` (500 obs) while primary model uses `dataListAgg` (250 obs).

```r
fit.lmmPP.half = lmerTest::lmer(PPvec ~ mapStimLevel + chanInCond +
    overallBlockType*pre_post*halfBlock +
    (1 + stim_linpoly|subjectNum:chanVec) + (1|subjectNum) + (1|subjectNum:blockVec),
    data=dataListAggHalf)
```

No significant effects involving halfBlock (three-way interaction F(3,417)=1.14, p=0.334; halfBlock main F(1,417)=0.01, p=0.909). EP amplitude is stable within blocks.

### Auxiliary Models

**`fit.lmmdiff`**: `absDiff ~ mapStimLevel + chanInCond + blockType + (1|subjectNum/chanVec)` -- models the absolute difference from baseline.

QQ plots for all models are saved as `pairedPulse_qq_*_allSubjs_AUC.png`.

### Workspace

All model objects and data are saved to `R_output/main_analysis_workspace.RData` for quick loading without rerunning the full script: `load(here("DBS_EP_PairedPulse","R_output","main_analysis_workspace.RData"))`.

### Plots

All EMM-related plots are saved as both PNG (300 dpi) and EPS in square dimensions (7x7 in).

1. Interaction plot: overallBlockType x pre_post with 95% CIs (dodged)
2. EMM by condition with CIs
3. Pre-post contrast within each condition (95% CI)
4. Standardized effect sizes (Cohen's d, total SD) for pre-post contrasts
5. Fixed-effect comparison across all model approaches
6. Cell-median percent difference from matched baseline (`across_subj_median_log_difference_model_data`) — uses the same cell medians as the statistical model, with percent difference computed from matched baseline cell medians

---

## 2. Single-Subject Anesthesia Analysis: `baseline_variability_3d413.R`

### Design

Single subject (3d413), 2 channels (4, 6), 5 blocks per channel (3 asleep, 2 awake), 4 stim levels (all levels included, `min_stim_level = 1`). 474 total trial-level observations (after 5-trial averaging, 10-1000 uV amplitude filter). Channels are recorded on alternating blocks (chan 6 → odd blocks, chan 4 → even blocks) due to alternating stimulation configurations.

**Block structure (channels alternate by block):**

| Channel | Block | blockType | Trials per stim level |
|---------|-------|-----------|-----------------------|
| 6 | 1 | asleep | ~12 |
| 4 | 2 | asleep | ~12 |
| 6 | 3 | awake | ~11-12 |
| 4 | 4 | awake | ~12 |
| 6 | 5 | awake | ~11-12 |
| 4 | 6 | awake | ~12 |
| 6 | 7 | asleep | ~12 |
| 4 | 8 | asleep | ~12 |
| 6 | 9 | asleep | ~12 |
| 4 | 10 | asleep | ~12 |

Per channel: 2 awake blocks, 3 asleep blocks, 5 blocks total.

### Models

Two LMMs are fit on trial-level data with block RE + uncorrelated random linear dose slope per block:

```r
fit.lmm1 = lmerTest::lmer(PPvec ~ mapStimLevel + blockType + chanVec +
    (1|blockVec) + (0+stim_linpoly|blockVec), data=dataList)  # additive
fit.lmm2 = lmerTest::lmer(PPvec ~ mapStimLevel * blockType + chanVec +
    (1|blockVec) + (0+stim_linpoly|blockVec), data=dataList)  # interaction
```

- `(1|blockVec)` -- random intercept per block (10 groups), captures overall amplitude differences between recording sessions
- `(0+stim_linpoly|blockVec)` -- random linear dose slope per block (uncorrelated with intercept), captures session-to-session variability in dose-response steepness. Fixes stim-level pseudoreplication (df drops from ~459 to ~23). The correlated version `(1+stim_linpoly|blockVec)` hit a boundary (correlation=1.00) with only 10 groups, so uncorrelated slopes are used.
- `chanVec` as fixed effect (only 2 channels, too few for RE)
- `blockType` (awake/asleep) is the effect of interest -- a between-block variable
- `mapStimLevel` as ordered factor (polynomial contrasts: .L, .Q, .C) -- a within-block variable
- `stim_linpoly` is the linear polynomial contrast extracted from `contr.poly(n_stim_levels)[, 1]`
- `trim_data = FALSE`, `log_data = FALSE`, amplitude filter 10-1000 uV

Log transformation was tested but produces catastrophic residual non-normality for this subject (skewness=-5.7, kurtosis=65) due to near-zero EP values creating a long left tail in log space. Raw-scale residuals are much better behaved. Log is needed for the multi-subject analysis (variance stabilization across subjects) but not for single-subject data with a more homogeneous variance structure.

### Degrees of Freedom and Satterthwaite Approximation

Satterthwaite df correctly reflect the hierarchy of the design:

| Effect | Type | DenDF | Interpretation |
|--------|------|-------|---------------|
| mapStimLevel (.L) | within-block (random slope) | ~9 | 10 blocks with random slope |
| mapStimLevel (.Q, .C) | within-block | ~452 | Residual-level (not absorbed by linear slope) |
| blockType | between-block | ~7 | 10 blocks - 3 between-block params (intercept + blockType + chanVec) |
| chanVec | between-block | ~7 | Same -- channels are assigned to alternating blocks |

### LMM Results

ANOVA (Type III, Satterthwaite):

| Effect | F | NumDF | DenDF | p |
|--------|---|-------|-------|---|
| mapStimLevel | 77.8 | 3 | 23.4 | 2.5e-12 |
| blockType | 0.09 | 1 | 7.0 | 0.774 |
| chanVec | 50.8 | 1 | 7.0 | < 0.001 |

Random effects: block intercept SD = 14.3, block slope SD = 44.1, residual SD = 41.3.

**No significant awake/asleep effect** (p = 0.774). Stimulation level drives EP amplitude (massive effect), and channel 6 has substantially larger EPs than channel 4 (p < 0.001).

### Permutation Tests (Trial-Level, Stratified)

Trial-level permutation tests complement the LMM by avoiding Satterthwaite assumptions. Tests are run **per channel** (channels are independent experiments) with 10,000 permutations.

**Approach:** Within each (channel, stim_level) stratum, awake/asleep labels are shuffled across trials. The stratified pooled test independently permutes within each stim level, then averages the per-stim-level median differences. This preserves stim-level balance in every permutation (important because stim level has a massive effect on amplitude).

**Trial-level exchangeability caveat:** The test shuffles individual trials, not blocks. With ICC = 0.074, trials within a block are slightly correlated, making the test mildly liberal. The effective sample size inflation factor is approximately `1 + (k-1) * ICC` where k ~ 12 trials/block, giving ~1.77x. This could make a nominal p = 0.05 actually ~0.07-0.08. Block-level permutation would be exact but has insufficient resolution: C(5,2) = 10 arrangements per channel, minimum p = 0.10.

**Per channel, per stim level:**

| Channel | Stim Level | Observed (uV) | p |
|---------|-----------|---------------|------|
| 6 | 1 | -13.4 | 0.111 |
| 6 | 2 | 11.9 | 0.537 |
| 6 | 3 | 14.3 | 0.385 |
| 6 | 4 | 18.6 | 0.279 |
| 4 | 1 | 3.3 | 0.532 |
| 4 | 2 | -0.8 | 0.946 |
| 4 | 3 | 12.7 | 0.262 |
| 4 | 4 | 5.9 | 0.303 |

**Stratified pooled (per channel):**

| Channel | Observed (uV) | p |
|---------|---------------|------|
| 6 | 7.9 | 0.30 |
| 4 | 5.3 | 0.17 |

All p-values are two-sided: `mean(abs(perm_diffs) >= abs(obs_diff))`.

**Consistent with LMM:** No significant awake/asleep effect in either channel.

### Effect Size Analysis

Cohen's d (awake vs asleep) computed per (channel, stim_level) on trial-level data. Raw uV scale (not log-transformed).

---

## 3. Single-Subject Conditioning Length Analysis: `length_conditioning_a23ed.R`

### Design

Single subject (a23ed), 1 channel (5), 4 blocks with **separate baselines**, 3 stim levels (after `min_stim_level = 2` filter). 134 total trial-level observations.
- Comparison 1: block 2 (baseline) vs block 3 (post 5-min A/B 200 ms)
- Comparison 2: block 5 (baseline) vs block 6 (post 15-min A/B 200 ms)

The two baselines are different measurement blocks separated in time. `trim_data = FALSE`, `log_data = FALSE` (raw uV for interpretability).

**Trial counts per cell:**

| Comparison | Stim Level | Pre (baseline) | Post (conditioning) | Total |
|-----------|-----------|----------------|---------------------|-------|
| 5-min | 2 | 10 | 11 | 21 |
| 5-min | 3 | 10 | 11 | 21 |
| 5-min | 4 | 9 | 11 | 20 |
| 15-min | 2 | 12 | 12 | 24 |
| 15-min | 3 | 12 | 12 | 24 |
| 15-min | 4 | 12 | 12 | 24 |

### Models

The mixed model approach is fundamentally limited here: with only 4 blocks, there are too few groups for a random effect (need >= 5-6). Block cannot be a fixed effect because it is perfectly confounded with `overallBlockType*pre_post`. The mixed models in the script (full interaction and reduced pre_post-only) are included for completeness but have convergence warnings and should be interpreted with caution.

Side-by-side comparison saved to `R_output/tab_model_a23ed_comparison.html`.

### Permutation Tests (Primary Statistical Inference)

Non-parametric permutation tests (10,000 shuffles) provide the primary inference. The test statistic is the **difference of medians**. All tests use trial-level permutation (shuffling block labels across trials). With only 2 blocks per comparison, block-level permutation is impossible (C(2,1) = 2 arrangements, minimum p = 1.0), so trial-level is the only feasible approach.

**Per-stim-level tests (Tests 1-3):** Permutation is within each stim level, preserving stim-level balance.

**Test 1 -- 5-min post vs baseline (per stim level):** Shuffles trials between blocks 2 and 3.

**Test 2 -- 15-min post vs baseline (per stim level):** Shuffles trials between blocks 5 and 6.

**Test 3 -- Per-stim-level interaction (15-min effect minus 5-min effect):** Independently shuffles within each comparison at each stim level, then computes the difference of differences.

Results (median-based, untrimmed):
- 5-min conditioning: marginal at stim 2 (p=0.055), not significant at stim 3 (p=0.75), significant at stim 4 (p=0.045)
- 15-min conditioning: highly significant at all stim levels (p<0.001, p<0.001, p=0.001)
- 15-min minus 5-min interaction: marginal at stim 2 (p=0.053), not significant at stim 3/4 (p=0.14, p=0.18)

**Stratified pooled tests:** These tests pool across stim levels while preserving stim-level balance. Within each permutation iteration, labels are shuffled independently within each stim level, per-stim-level median differences are computed, then averaged. This avoids the problem with naive pooling (mixing trials across stim levels), where random stim-level imbalance between shuffled groups creates spurious differences that widen the null distribution and make the test conservative.

**Stratified pooled interaction (15-min effect minus 5-min effect):**
- Observed: 60.7 uV, **p = 0.0034**

**Stratified pooled per-condition:**

| Condition | Observed (uV) | p |
|-----------|---------------|------|
| A/B 200 ms 5 minutes | 16.5 | 0.013 |
| A/B 200 ms 15 minutes | 77.2 | < 0.0001 |

Null distribution histograms are saved for the pooled interaction and per-condition tests.

P-values are two-sided: `mean(abs(perm_diffs) >= abs(obs_diff))`.

**Important caveats:**

1. **Trial exchangeability.** The permutation test shuffles individual trials between blocks, assuming they are exchangeable under the null. This is the same independence assumption as a t-test or lm -- if trials within a block are temporally correlated, the permutation null distribution is too narrow and p-values may be anti-conservative. The 5-trial averaging upstream partially mitigates temporal autocorrelation but does not eliminate it. With only 2 blocks per comparison, there is no way to estimate within-block ICC for this subject.

2. **Permutation space.** With ~10 trials per block per stim level, the number of unique permutations per cell is C(~20, ~10) ~ 184,756. The minimum achievable p-value is ~5e-06. The reported p<0.001 values are well above this floor. The stratified pooled tests combine 3 independent stim-level shuffles, giving an even richer permutation space (product of per-stim-level spaces).

3. **Stratified vs non-stratified pooling.** Per-stim-level tests and stratified pooled tests never mix trials from different stimulation intensities -- the stim-level balance is preserved by construction. This is important because stim level has a large effect on EP amplitude.

### Effect Size Analysis

Cohen's d (post vs baseline) computed per (stim_level, conditioning_length) on trial-level data. PPvec is NOT log-transformed (`log_data=FALSE`), so d is on original uV scale. Note this makes effect sizes not directly comparable to the main multi-subject analysis (which uses log scale).

### Plots

1. Per-subject scatter plots by block
2. Across-condition scatter and summary plots
3. EMM interaction plot (overallBlockType x pre_post)
4. Cohen's d bar plot with CI error bars and small/medium/large reference lines
5. Stratified permutation null distribution histogram: 5-min post vs baseline
6. Stratified permutation null distribution histogram: 15-min post vs baseline
7. Stratified permutation null distribution histogram: 15-min effect minus 5-min effect (interaction)
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
| No permutation tests for 3d413 | Medium | **Fixed** -- trial-level stratified permutation tests per channel |
| a23ed pooled permutation not stratified by stim level | Medium | **Fixed** -- stratified pooling (independent permutation within each stim level) |
| Trimming + median redundancy | Medium | **Fixed** -- `trim_data = FALSE` with median |
| Effect size using residual-only SD | Medium | **Fixed** -- total SD for comparability |
| CIs using single approximate t_crit | Low-Med | **Fixed** -- emmeans-native CIs per contrast |
| AR(1) degenerate/overparameterized | Medium | **Fixed** -- switched to lmer |
| `emm_options` global side effect | Medium | **Fixed** -- scoped per model |
| `as_data_frame()` deprecation | Low | Noted -- should migrate to `as_tibble()` |

## Known Limitations

- **Interaction is not significant.** The conditioning x pre_post interaction (p=0.251) does not reach significance. However, the conditioning protocol main effect is now significant (p=0.044), driven by A/B 200 vs A/B 25 (Tukey p=0.030). It was previously marginal (p=0.068 with min_stim_level=3, no amplitude filter) but is sensitive to the number of conditions included, aggregation method, RE structure, and data filtering.
- **a23ed: 4 blocks is insufficient for mixed models.** Block RE requires >= 5-6 groups. All mixed models for a23ed have convergence warnings. Permutation tests are the primary analysis, with the caveat that they assume trial exchangeability.
- **Trial-level permutation assumes exchangeability.** Both single-subject analyses (3d413, a23ed) use trial-level permutation -- shuffling individual trial labels rather than block labels. This assumes trials are exchangeable across blocks under the null, which is violated when there are block effects beyond the condition of interest (ICC ~ 0.07 for 3d413; unestimable for a23ed with only 2 blocks per comparison). The effective sample size inflation factor is approximately `1 + (k-1) * ICC` where k is the cluster size. With ICC = 0.07 and k ~ 12, this is ~1.77x, potentially making a nominal p = 0.05 actually ~0.07-0.08. Block-level permutation would be exact but is infeasible: C(5,2) = 10 per channel for 3d413 (minimum p = 0.10), C(2,1) = 2 for a23ed (minimum p = 1.0). Trial-level permutation is the only approach with adequate resolution.
- **3d413: Low power for awake/asleep effect.** With only ~7 Satterthwaite df for blockType (10 blocks - 3 between-block params), and 4 awake vs 6 asleep blocks, the LMM has limited power. The permutation tests confirm the non-significant result but share the same fundamental sample size limitation.
- **3d413: Satterthwaite approximation with few groups.** Satterthwaite df depend on well-estimated variance components. With 10 blocks, the block variance has high uncertainty, and normality of block random effects is unverifiable. The approximation may be mildly anti-conservative. Permutation tests sidestep these assumptions entirely.
- **Cohen's d on different scales.** Main script uses log(uV), 3d413 and a23ed use raw uV. Effect sizes are not directly comparable between analyses.
- **Small sample size.** 13 subjects, 51 recording sessions, 250 cell medians. Power for block-level effects is limited by the ~51 block groups.
- **Disease and Brodmann area underpowered.** 4 MD, 9 PD patients. baLabel has 4 levels (3 Brodmann areas + "Unknown" for a23ed) across 17 channels; tested and not significant (F(3,11)=0.20, p=0.90).
- **Unequal 5-trial averaging.** The last averaged trial in a block may be the average of 1-4 trials (when trial count is not divisible by 5), giving it higher variance. Cell-median aggregation partially mitigates this.
- **mapStimLevel as ordered factor.** Polynomial contrasts assume equally spaced levels, but actual stimulation currents (mA) are not equally spaced across subjects. The contrasts are meaningful for rank order only, not physical units.

## Contrast Coding

R's default contrast settings apply throughout (no scripts override `options(contrasts=...)`):
- **Unordered factors** (`as.factor()`) → `contr.treatment` (dummy coding, compare each level to reference)
- **Ordered factors** (`as.ordered()`) → `contr.poly` (orthogonal polynomial: .L linear, .Q quadratic, .C cubic, etc.)

Contrasts are per-variable properties in R — setting one variable as ordered does not affect other variables. No global reset is needed.

| Variable | Script(s) | Type | Contrasts | Levels | Notes |
|----------|-----------|------|-----------|--------|-------|
| `mapStimLevel` | all three | `as.ordered()` | polynomial | up to 4 (all levels included, `min_stim_level = 1`) | Ordinal dose-response. .L tests linear trend, .Q tests curvature, .C tests cubic. Polynomial contrasts assume equally-spaced levels; mapped integers are equally spaced but underlying mA values differ per subject. |
| `blockType` | 3d413, a23ed | `as.factor()` | treatment | 2 | awake/asleep or baseline/condition. With 2 levels, treatment and polynomial contrasts give identical F-tests (1 df). |
| `chanVec` | 3d413 | `as.factor()` | treatment | 2 | Channels 4, 6. 1 df — contrast type irrelevant for F-test. |
| `overallBlockType` | main, a23ed | `as.factor()` | treatment | 4 (main), 2 (a23ed) | Conditioning protocol. Treatment contrasts compare each to reference; doesn't affect Type III F-test (invariant to parameterization). Pairwise comparisons via `emmeans` are contrast-agnostic. |
| `pre_post` | main, a23ed | character (auto-converted by lmer) | treatment | 2 | 1 df — contrast type irrelevant. |
| `chanInCond` | main | `as.factor()` | treatment | 2 | Whether channel was in conditioning pair. 1 df. |
| `subjectNum` | main | `as.factor()` | treatment | 9 | Only appears in RE terms `(1\|subjectNum/...)`, not as fixed effect — contrasts are unused. |
| `blockVec` | all three | `as.factor()` | treatment | varies | Only appears in RE terms `(1\|blockVec)` — contrasts are unused. |

**Key points:**
1. **Type III ANOVA is invariant to contrast coding.** The F-test for each fixed effect tests whether *any* level means differ, regardless of parameterization. Contrast coding only affects individual coefficient estimates.
2. **All 2-level factors are contrast-agnostic.** With 1 df, treatment, sum, and polynomial contrasts produce identical F-tests and p-values.
3. **`mapStimLevel` polynomial contrasts are scientifically appropriate** for an ordinal dose-response variable. The .L (linear) component captures the primary effect of interest (amplitude increases with stimulation intensity). The .Q (quadratic) component tests for diminishing returns or saturation.
4. **No manual contrast manipulation.** No script calls `contrasts(var) <- ...`. All contrasts come from R's default behavior based on `as.factor()` vs `as.ordered()`.

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
- `subj_a23ed_perm_interaction.png/eps` -- stratified permutation null: 15-min vs 5-min interaction
- `subj_a23ed_perm_a_b_200_ms_5_minutes.png/eps` -- stratified permutation null: 5-min conditioning
- `subj_a23ed_perm_a_b_200_ms_15_minutes.png/eps` -- stratified permutation null: 15-min conditioning
- `subj_a23ed_two_panel_prepost.png/eps` -- two-panel pre/post with separate baselines
- `subj_a23ed_effect_size_conditioning_length.png/eps` -- Cohen's d bar plot
- `subj_3d413_perm_awake_asleep_chan*.png/eps` -- stratified permutation null: awake vs asleep per channel
- `subj_3d413_effect_size_awake_asleep.png/eps` -- Cohen's d bar plot (raw uV scale)
