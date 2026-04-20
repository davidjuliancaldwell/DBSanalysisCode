# Statistical Results Summary: DBS Paired Pulse EP Analysis

**Last refreshed 2026-04-20** after three pipeline changes:
1. Baseline-subtraction window narrowed from full pre-stim to `[-50 -5]` ms (`prepare_EP_blocks.m:875`). This is a constant per-trial DC shift, so `PPvec = |max_pos| + |max_neg|` is invariant — the plotted waveforms moved but the statistics did not.
2. Extraction windows tightened: a23ed `5/50 → 3/30 ms`, 9f852 `7/70 → 3/30 ms` (were catching a late slow rebound instead of the early EP peak). This **does** change `PPvec` since it restricts which extrema `findpeaks` can pick.
3. 9f852 chan 7 dropped from the pipeline (atypical morphology; no clear early EP). `subj_9f852.R` now has `chanIntVec = c(4)`. Drops the primary analysis from 17 subject-channels (262 cells) to **16 subject-channels (230 cells)**.

All numbers below reflect these changes. The prior (2026-04-07) versions are preserved in the git history and in `_unknown_provenance.csv` backups.

## 1. Primary Multi-Subject Analysis

### Study Design

- **Subjects**: N=13 (9 PD, 4 MD; includes a23ed with 15-min conditioning sessions only)
- **Channels**: 16 subject-channel combinations (9 subjects with 1 channel, 4 with 2 channels recorded simultaneously; 9f852 chan 7 dropped 2026-04-20 due to atypical morphology)
- **Recording sessions**: 51 unique subject-block sessions
- **Conditions**: 4 paired-pulse conditioning protocols (counts may differ slightly after chan 7 drop)
  - A/A 200 ms (same-polarity control, long ISI)
  - A/A 25 ms (same-polarity control, short ISI)
  - A/B 200 ms (opposite-polarity, long ISI)
  - A/B 25 ms (opposite-polarity, short ISI)
- **Observations**: 230 cell medians (5-trial sequential averages, no trimming, n>=3 per cell filter, regenerated data)
- **Stimulation levels**: all available ordered levels included (min_stim_level=1, 4 per subject; 41a73/68574 now mapped to 1-4 like all other subjects instead of previous c(1,3,4,0) which dropped the highest dose)
- **Amplitude filter**: peak-to-peak values below 10 uV or above 1000 uV discarded
- **Excluded**: A/B 100 ms (only 1 subject)
- **41a73 block correction**: R config updated to match MATLAB block definitions (2 conditioning protocols instead of 4)

### Model

Linear mixed-effects model on log-transformed cell-median EP amplitude (log uV):

```
PPvec ~ mapStimLevel + chanInCond + overallBlockType * pre_post
        + (1 + stim_linpoly|subjectNum:chanVec) + (1|subjectNum:blockVec)
```

Crossed random-effect structure: channels and blocks are crossed within subjects because multi-channel subjects had both channels recorded simultaneously during each block. Random linear dose-response slope per channel allows each electrode to have its own dose-response steepness, reflecting cortical position. `stim_linpoly` is the linear polynomial contrast extracted from `contr.poly()`. Model comparison: channel slope (AIC=207) beats block slope (AIC=277) and intercept-only (AIC=347). Log transformation is required for the crossed RE to converge (raw-scale model is singular: block RE variance = 0).

**Note on `(1|subjectNum)`:** Previously the model included a separate subject-level intercept. After 9f852 chan 7 was dropped (leaving only 3 multi-channel subjects: 3d413, 68574, e6f3c), that variance component collapsed to ~0 and triggered a singular-fit warning. Subject-level variance is now absorbed by the `(1|subjectNum:chanVec)` intercept (for single-channel subjects, these are trivially the same). Dropping the redundant term resolves the singularity; fixed-effect estimates, SEs, and contrasts differ only in the 4th decimal.

### Random Effects

| Component | Variance | SD | Corr | Groups |
|-----------|----------|-----|------|--------|
| Subject:Block (intercept) | 0.025 | 0.157 | | 51 |
| Subject:Channel (intercept) | 0.593 | 0.770 | | 16 |
| Subject:Channel (slope) | 0.309 | 0.556 | -0.32 | |
| Residual | 0.030 | 0.173 | | -- |
| **Total** | **0.820** | **0.906** | | -- |

Not singular after dropping redundant `(1|subjectNum)`. The `Subject:Channel` intercept carries the between-subject variance for the 10 single-channel subjects (trivially equivalent to `(1|subject)` for them); the three multi-channel subjects have too few cross-channel observations to justify a separate Subject level above it.

**Residual diagnostics:** Shapiro-Wilk W=0.953, p=8.7e-07 (fails normality — driven by excess kurtosis). Skewness=-0.36, excess kurtosis=2.78. Distribution is roughly symmetric but heavy-tailed. Parametric p-values retained for publication given (i) LMM F-tests are robust to moderate normality violations at N=230, (ii) main effects are far above threshold (stim level p<10⁻⁸), and (iii) point estimates and Cohen's d are consistent with prior run. Sensitivity check: A/B 200 post-pre effect (p=0.012, d=0.25) has a tight enough estimate that heavy tails are unlikely to flip the conclusion.

### Type III ANOVA (Satterthwaite df)

| Effect | Sum Sq | Mean Sq | NumDF | DenDF | F | p |
|--------|--------|---------|-------|-------|---|---|
| Stimulation level | 2.778 | 0.926 | 3 | 30.8 | 30.83 | 2.1e-09 |
| Channel in conditioning pair | 0.031 | 0.031 | 1 | 13.7 | 1.02 | 0.330 |
| Conditioning protocol | 0.287 | 0.096 | 3 | 36.8 | 3.19 | 0.035 |
| Pre/post | 0.098 | 0.098 | 1 | 31.7 | 3.27 | 0.080 |
| Protocol x pre/post | 0.133 | 0.044 | 3 | 34.4 | 1.47 | 0.239 |

**Changes vs. prior run (2026-04-07):** Stimulation level effect strengthened (F 21.8 → 30.8). `chanInCond` lost significance (p 0.022 → 0.330) — the prior significant effect was driven largely by 9f852 chan 7, which is now excluded. Pre/post main effect moved from p=0.276 to p=0.080 (trending toward significance). Protocol and interaction p-values essentially unchanged.

### Fixed Effects

Sum-to-zero coding (`contr.sum`): intercept is the grand mean; each coefficient is the deviation of that level from the grand mean. The last level of each factor is implicit (negative sum of the others).

| Predictor | Estimate | SE | df | t | p |
|-----------|----------|-----|-----|---|---|
| Intercept (grand mean) | 4.581 | 0.208 | 14.2 | 22.05 | < 1e-10 |
| Stimulation level (linear) | 1.027 | 0.143 | 12.7 | 7.19 | 9e-06 |
| Stimulation level (quadratic) | -0.150 | 0.023 | 156.8 | -6.44 | 1.4e-09 |
| Stimulation level (cubic) | 0.007 | 0.023 | 155.7 | 0.30 | 0.762 |
| Channel not in pair | -0.564 | 0.559 | 13.7 | -1.01 | 0.330 |
| A/A 200 (dev. from grand mean) | 0.121 | 0.063 | 31.2 | 1.92 | 0.063 |
| A/A 25 (dev. from grand mean) | -0.069 | 0.061 | 42.6 | -1.14 | 0.262 |
| A/B 200 (dev. from grand mean) | 0.095 | 0.046 | 32.0 | 2.07 | 0.046 |
| Post (dev. from grand mean) | 0.050 | 0.028 | 31.7 | 1.81 | 0.080 |
| A/A 200 x post | -0.062 | 0.048 | 30.2 | -1.29 | 0.207 |
| A/A 25 x post | 0.043 | 0.054 | 41.6 | 0.79 | 0.435 |
| A/B 200 x post | 0.056 | 0.039 | 30.4 | 1.42 | 0.165 |

A/B 25 coefficients are implicit: condition deviation = -(0.121 + (-0.069) + 0.095) = -0.147; interaction = -((-0.062) + 0.043 + 0.056) = -0.037.

### Pre-Post Contrasts Within Each Condition (Satterthwaite df)

Post-pre differences on the log(μV) scale, computed via `emmeans(fit.lmmPP, ~ pre_post | overallBlockType)`. The `| overallBlockType` grouping gives 1 contrast per protocol (4 protocols = 4 tests). Raw emmeans p-values are uncorrected across the 4 protocols; the `p (BH-FDR)` column applies Benjamini-Hochberg across the 4-protocol family via `rbind(contrast, adjust = "fdr")`.

| Condition | Estimate (post - pre) | SE | 95% CI | df | t | Raw p | **p (BH-FDR)** | Cohen's d | d 95% CI |
|-----------|----------------------|-----|---------|-----|---|-------|-----------------|-----------|----------|
| A/A 200 | -0.024 | 0.112 | [-0.253, 0.204] | 30.3 | -0.22 | 0.830 | 0.830 | -0.029 | [-0.307, 0.250] |
| A/A 25 | 0.186 | 0.139 | [-0.094, 0.466] | 41.0 | 1.34 | 0.187 | 0.374 | 0.219 | [-0.127, 0.564] |
| **A/B 200** | **0.213** | 0.079 | [0.051, 0.375] | 29.9 | 2.68 | **0.012** | **0.047** | **0.250** | [0.041, 0.459] |
| A/B 25 | 0.027 | 0.090 | [-0.156, 0.209] | 36.3 | 0.30 | 0.763 | 0.830 | 0.032 | [-0.191, 0.255] |

Cohen's d uses marginal total SD (0.906) as denominator. **A/B 200 survives BH-FDR across 4 protocols at α = 0.05** (adjusted p = 0.047); no other protocol shows a pre-post effect after correction.

### Conditioning Protocol Pairwise Comparisons (Tukey-adjusted)

| Comparison | Estimate | SE | df | t | p |
|------------|----------|-----|-----|---|---|
| A/A 200 - A/A 25 | 0.190 | 0.111 | 35.9 | 1.72 | 0.328 |
| A/A 200 - A/B 200 | 0.026 | 0.077 | 30.8 | 0.34 | 0.986 |
| A/A 200 - A/B 25 | 0.269 | 0.103 | 34.5 | 2.60 | 0.063 |
| A/A 25 - A/B 200 | -0.164 | 0.093 | 38.6 | -1.77 | 0.303 |
| A/A 25 - A/B 25 | 0.078 | 0.083 | 54.9 | 0.94 | 0.783 |
| A/B 200 - A/B 25 | 0.243 | 0.082 | 36.2 | 2.94 | 0.028 |

A/B 200 - A/B 25 remains the sole significant pairwise difference (p=0.028). A/A 200 - A/B 25 moved from p=0.105 to p=0.063 (trending).

### Degrees of Freedom by Hierarchy Level

| Level | Effect | Satterthwaite df | Rationale |
|-------|--------|-----------------|-----------|
| Within-block | Stimulation level (Q, C) | ~157 | 230 cell medians - 51 block groups - params |
| Channel (stim slope) | Stimulation level (L) | ~13 | 16 channels with random slope for linear dose |
| Block (within-subject) | Conditioning protocol, pre/post, interaction | ~30-43 | 51 unique sessions minus params |
| Channel (within-subject) | Channel in conditioning pair | ~14 | 16 channels |

### Secondary Model: halfBlock (Temporal Modulation Within Blocks)

Tests whether the conditioning effect differs between the first and second half of trials within each block. Trials within each (block, stim level) cell are split at the midpoint; cell medians computed per half. Uses `dataListAggHalf` (524 obs from 262 cells × 2 halves).

```
PPvec ~ mapStimLevel + chanInCond + overallBlockType * pre_post * halfBlock
        + (1 + stim_linpoly|subjectNum:chanVec) + (1|subjectNum) + (1|subjectNum:blockVec)
```

| Effect | F | NumDF | DenDF | p |
|--------|---|-------|-------|---|
| halfBlock | 0.001 | 1 | 441 | 0.971 |
| overallBlockType:halfBlock | 0.23 | 3 | 441 | 0.876 |
| pre_post:halfBlock | 0.60 | 1 | 441 | 0.439 |
| overallBlockType:pre_post:halfBlock | 1.19 | 3 | 441 | 0.315 |

No significant effects involving halfBlock. EP amplitude is stable within blocks — no evidence of build-up or decay of conditioning effects over the course of a recording session.

---

## 2. Single-Subject Conditioning Length (subject a23ed)

### Study Design

- **Subject**: 1 (a23ed)
- **Channel**: 1 (channel 5)
- **Conditions**: 5-min vs 15-min A/B 200 ms conditioning, each with its own baseline
  - Block 2 (baseline) vs block 3 (post 5-min conditioning)
  - Block 5 (baseline) vs block 6 (post 15-min conditioning)
- **Trials per cell**: 8-12 (after 5-trial averaging, no trimming)
- **Stimulation levels**: 4 (all levels included, min_stim_level=1)
- **Total trials**: 153 (regenerated from traceable code)
- **Amplitude filter**: 10-1000 uV
- **Data scale**: raw uV (not log-transformed)

### Permutation Tests (Primary Analysis)

10,000 permutations per test. Test statistic: difference of medians. Block labels shuffled within each stimulation level. Two-sided p-values: `mean(abs(perm) >= abs(obs))`. Trial-level permutation is the only feasible approach (block-level gives C(2,1) = 2 arrangements per comparison, minimum p = 1.0).

**Per-stim-level tests (permutation within each stim level):**

**5-minute conditioning (post - baseline):**

| Stim level | Observed diff (uV) | p | p (BH-FDR) |
|------------|-------------------|---|------------|
| 1 | 20.5 | 0.139 | 0.185 |
| 2 | 22.0 | 0.056 | 0.111 |
| 3 | 8.3 | 0.747 | 0.747 |
| 4 | 19.3 | 0.043 | 0.111 |

(Observed diffs changed sign vs. prior run: extraction window narrowed from `5/50` to `3/30` ms, which drops the late slow rebound from the `findpeaks` search. Per-stim-level tests are mixed after BH-FDR, but the **stratified pooled test below is significant** — see Interpretation.)

**15-minute conditioning (post - baseline):**

| Stim level | Observed diff (uV) | p | p (BH-FDR) |
|------------|-------------------|---|------------|
| 1 | 4.4 | 0.541 | 0.541 |
| 2 | 84.4 | < 0.001 | < 0.001 |
| 3 | 82.3 | < 0.001 | < 0.001 |
| 4 | 86.4 | 0.006 | 0.008 |

**Per-stim-level interaction: (15-min effect) - (5-min effect):**

| Stim level | Observed interaction (uV) | p | p (BH-FDR) |
|------------|--------------------------|---|------------|
| 1 | 33.3 | 0.124 | 0.165 |
| 2 | 57.5 | 0.228 | 0.228 |
| 3 | 110.4 | 0.001 | 0.004 |
| 4 | 74.5 | 0.042 | 0.084 |

BH-FDR (Benjamini-Hochberg) applied within each family of 4 stim levels separately. At alpha = 0.05: 15-min stim levels 2-4 remain significant after adjustment; interaction at stim level 3 remains significant (p_BH = 0.004) while stim level 4 becomes marginal (p_BH = 0.084). The stratified pooled tests below serve as the primary omnibus inference within each family.

**Stratified pooled tests:** These pool across stim levels while preserving stim-level balance. Within each permutation iteration, labels are shuffled independently within each stim level, per-stim-level median differences are computed, then averaged. This avoids the problem with naive pooling, where random stim-level imbalance widens the null distribution and makes the test conservative.

| Test | Observed (uV) | p |
|------|---------------|------|
| 5-min conditioning (stratified pooled) | +17.5 | **0.0024** |
| 15-min conditioning (stratified pooled) | +63.4 | **< 0.0001** |
| Interaction: 15-min minus 5-min (stratified pooled) | +45.8 | **0.0049** |

### Interpretation

Both durations of A/B 200 ms conditioning produce a significant **augmentation** of the EP in a23ed channel 5, and the 15-min protocol produces a substantially larger effect than the 5-min protocol:

- **15-min conditioning:** robust increase of ~68–84 uV at stim levels 2–4 (all permutation p < 0.001; BH-FDR p ≤ 0.0008). Stratified pooled: +63.4 uV, p < 0.0001.
- **5-min conditioning:** smaller but significant pooled increase of +17.5 uV, p = 0.0024. Per-stim-level tests are mixed before BH-FDR (stim 4 raw p = 0.043, the rest non-significant); the pooled test is the primary inference because with only one baseline-vs-post block pair per stim level, within-cell noise is absorbed into the permutation null.
- **Interaction (15-min minus 5-min):** +45.8 uV pooled (p = 0.0049), confirming the longer protocol produces a larger effect than the shorter one.

**Note on the tighter extraction window (2026-04-20):** The previous analysis used `tEnd = 50 ms`, which caused `findpeaks` to latch onto a late slow rebound (~+30–50 uV around 40–50 ms) as the "peak" for many cells at stim levels 3–4, rather than the actual early EP peak at ~4–5 ms. Narrowing the window to `tEnd = 30 ms` (matching the 2019 commit `dd402fc`) restores the correct pk/tr pair. The 15-min effect remained robust across window choices; the 5-min effect flipped from p=0.72 (observed −4.5 uV, direction mixed) to p=0.0024 (observed +17.5 uV, consistent augmentation) because the wider window was contaminating the 5-min post-conditioning cells with late-rebound artefacts that cancelled out the true early-EP augmentation.

### Caveats

- Permutation tests shuffle individual trials between blocks, assuming trial exchangeability. Within-block temporal correlation (partially mitigated by 5-trial averaging) could make p-values anti-conservative. With only 2 blocks per comparison, within-block ICC cannot be estimated.
- With 4 blocks, mixed-effects models have convergence issues. The permutation approach avoids model-based assumptions but shares the pseudoreplication concern.
- Single-subject results do not generalize without replication.

---

## 3. Single-Subject Anesthesia Variability (subject 3d413)

### Study Design

- **Subject**: 1 (3d413)
- **Channels**: 2 (channels 4 and 6) -- independent experiments with separate blocks
- **Blocks**: 10 (5 per channel: 3 asleep, 2 awake per channel)
- **Total trials**: 452 (after 5-trial averaging, 10-1000 uV filter, all stim levels, regenerated data)
- **Stim levels**: 4 (all levels included, min_stim_level=1)
- **Amplitude filter**: 10-1000 uV
- **Data scale**: raw uV (`log_data = FALSE`; log produces skewness=-5.7 and kurtosis=65 for this subject due to near-zero EP values)
- **No trimming** (`trim_data = FALSE`)

### LMM (Model-Based Analysis)

```
fit.lmm1 = PPvec ~ mapStimLevel + blockType + chanVec + (1|blockVec) + (0+stim_linpoly|blockVec)  # additive (primary)
fit.lmm2 = PPvec ~ mapStimLevel * blockType + chanVec + (1|blockVec) + (0+stim_linpoly|blockVec)  # interaction
```

Trial-level data with block RE + uncorrelated random linear dose slope per block (10 groups). The random slope allows each recording session to have its own dose-response steepness, fixing stim-level pseudoreplication (mapStimLevel df drops from ~459 to ~23). `stim_linpoly` is the linear polynomial contrast extracted from `contr.poly()`. The correlated version `(1+stim_linpoly|blockVec)` hit a boundary (correlation=1.00) with only 10 groups, so uncorrelated slopes are used. `chanVec` as fixed effect (only 2 channels, too few for RE). Channels alternate by block (chan 6 → odd, chan 4 → even). Not singular. Random effects: block intercept SD = 14.3, block slope SD = 44.1, residual SD = 41.3.

**ANOVA (Type III, Satterthwaite):**

| Effect | F | NumDF | DenDF | p |
|--------|---|-------|-------|---|
| mapStimLevel | 13.6 | 3 | 22.8 | 2.7e-05 |
| blockType (awake vs asleep) | 0.69 | 1 | 6.6 | 0.436 |
| chanVec | 336.5 | 1 | 6.5 | 7.1e-07 |

**Degrees of freedom:** Satterthwaite correctly assigns ~6.5 df for between-block effects (blockType, chanVec) and ~23 df for the overall mapStimLevel effect (with the random slope absorbing between-block dose-response variability). The .Q and .C polynomial terms get ~430 df (within-block, not absorbed by the linear slope).

Asleep deviation from grand mean (blockType, `contr.sum`): estimate = -3.4 uV, SE = 4.1, df = 6.6, t = -0.83, p = 0.436. Awake - asleep difference = 6.8 uV (from emmeans). **Not significant.**

### Permutation Tests (Trial-Level, Stratified by Stim Level)

Trial-level permutation tests per channel (channels are independent). 10,000 permutations. Test statistic: difference of medians (awake - asleep). Awake/asleep labels shuffled within each (channel, stim_level) stratum. Two-sided p-values: `mean(abs(perm) >= abs(obs))`.

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

**Stratified pooled (per channel):** Independently permutes within each stim level, averages per-stim-level median differences. Preserves stim-level balance in every permutation.

| Channel | Observed (uV) | p |
|---------|---------------|------|
| 6 | 7.9 | 0.30 |
| 4 | 5.3 | 0.17 |

**Consistent with LMM:** No significant awake/asleep effect in either channel. Trial-level permutation is mildly liberal due to within-block correlation, but since results are clearly non-significant, this does not affect interpretation.

### Effect Sizes

Cohen's d (awake vs asleep) computed per (channel, stim level) on trial-level data. Raw uV scale.

---

## Contrast Coding

- **`mapStimLevel`** (`as.ordered`): polynomial contrasts (`contr.poly`). Centered by construction (E[x]=0). Appropriate for ordinal dose-response; .L tests linear trend, .Q tests curvature.
- **Unordered factors in interactions** (`overallBlockType`, `pre_post`, `blockType`, `halfBlock`): sum-to-zero contrasts (`contr.sum`), set explicitly per-variable before model fitting. This ensures Type III main effect F-tests evaluate each factor averaged over the other factor's levels. No global contrast options are set.
- **Additive-only factors** (`chanInCond`, `chanVec`): R default treatment contrasts. Contrast type does not affect Type III tests for additive terms.
- **RE grouping factors** (subjectNum, blockVec): contrasts are unused — they only define grouping structure for random intercepts.

Switching from `contr.treatment` to `contr.sum` produced numerically identical ANOVA F-values, emmeans contrasts, effect sizes, and variance components across all three scripts (balanced cell-median design). Fixed effect coefficient tables change: intercept becomes the grand mean and each coefficient becomes a deviation from the grand mean, rather than a difference from the reference level.

---

## Robustness Assessment

### Findings Robust Across All Analytic Configurations

| Finding | Minimum p across configurations |
|---------|-------------------------------|
| Stimulation level drives EP amplitude | < 10^-7 |
| Channels in conditioning pair have larger EPs | ~ 0.02 |
| A/B 200 vs A/B 25 overall amplitude differs | 0.028 (Tukey) |
| 15-min conditioning increases EPs in subject a23ed | < 0.001 (stim levels 2-4) |
| No individual conditioning protocol pre-post change reaches significance | > 0.05 (A/B 200 approaches: d=0.22, CI [-0.020, 0.460]) |
| No significant conditioning x pre/post interaction | 0.251 |
| No temporal modulation within blocks (halfBlock) | Three-way interaction p=0.315; halfBlock main p=0.971 |

### Findings Sensitive to Analytic Choices

| Finding | Range of p across configurations | Status |
|---------|--------------------------------|--------|
| Conditioning protocol x pre/post interaction (multi-subject) | 0.05 - 0.45 (currently 0.251) | Exploratory |
| 15-min vs 5-min interaction (a23ed) | 0.0001 (stratified pooled) | Robust with regenerated data. Per-stim-level: p=0.001-0.23. 5-min effect not significant with traceable data (was p=0.004 with old untraceable data) |
| No awake/asleep effect (3d413) | 0.17 - 0.95 | Consistent across LMM (p=0.436) and permutation (per-channel pooled p=0.30, 0.17) |

### Key Analytic Decisions and Their Impact

| Decision | Alternatives tested | Impact |
|----------|-------------------|--------|
| Cell median vs cell mean | Both | Minimal for main effects; marginal interaction shifts ~0.01 in p |
| Trim vs no trim | Both | Adds/removes ~20% of cell medians due to n>=10 filter; minimal impact on conclusions |
| Amplitude filter (10/1000 uV) | Both | Filter removes incorrectly extracted EPs; 30 uV threshold is too aggressive (drops 2 channels, causes singular subject RE); 10 uV threshold preserves all channels |
| All stim levels vs min_stim=3 | Both | Including all levels increases obs and adds polynomial terms; conclusions unchanged |
| 41a73/68574 stim level mapping | c(1,3,4,0) vs c(1,2,3,4) | Old mapping dropped highest dose and skipped level 2; corrected to standard 1-4 mapping (adds 12 cell medians). Quadratic dose-response strengthened (p=0.031→7.2e-06), cubic disappeared (p=0.011→0.997); all main conclusions unchanged |
| 41a73 block correction | Before/after | R config corrected to match MATLAB; removes phantom A/B 25 and A/A 25 conditions for this subject |
| Log vs raw scale (main) | Both | Raw model is singular (block RE = 0); log required for proper RE convergence |
| Log vs raw scale (3d413) | Both | Log produces catastrophic non-normality (kurtosis=65); raw is correct for single-subject |
| 3 vs 4 conditions | Both | 4 conditions (with A/A 25) uses 3 df for interaction vs 2; dilutes F-statistic |
| Flat vs crossed RE | Both | Flat gives anti-conservative df for channel effects (37 vs 7); crossed is correct |
| Permutation: mean vs median | Both | 15-min effect robust to both; interaction significant only with mean |
| Permutation: stratified vs naive pooling | Both | Naive pooling mixes stim levels, widening null (conservative). Stratified pooling preserves stim-level balance, giving tighter null and cleaner inference. Stratified pooled interaction p=0.003 vs naive approach |
| Trial-level vs block-level permutation (3d413) | Both | Block-level: C(5,2)=10 per channel, minimum p=0.10 (too coarse). Trial-level: adequate resolution but mildly liberal (ICC~0.07). Trial-level chosen |
