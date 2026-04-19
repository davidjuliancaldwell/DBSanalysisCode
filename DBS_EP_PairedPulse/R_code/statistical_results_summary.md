# Statistical Results Summary: DBS Paired Pulse EP Analysis

## 1. Primary Multi-Subject Analysis

### Study Design

- **Subjects**: N=13 (9 PD, 4 MD; includes a23ed with 15-min conditioning sessions only)
- **Channels**: 17 subject-channel combinations (9 subjects with 1 channel, 4 with 2 channels recorded simultaneously)
- **Recording sessions**: 51 unique subject-block sessions
- **Conditions**: 4 paired-pulse conditioning protocols
  - A/A 200 ms (same-polarity control, long ISI): 2 subjects, 4 blocks
  - A/A 25 ms (same-polarity control, short ISI): 3 subjects, 6 blocks
  - A/B 200 ms (opposite-polarity, long ISI): 8 subjects, 16 blocks
  - A/B 25 ms (opposite-polarity, short ISI): 7 subjects, 14 blocks
- **Observations**: 262 cell medians (5-trial sequential averages, no trimming, n>=3 per cell filter, regenerated data)
- **Stimulation levels**: all available ordered levels included (min_stim_level=1, 4 per subject; 41a73/68574 now mapped to 1-4 like all other subjects instead of previous c(1,3,4,0) which dropped the highest dose)
- **Amplitude filter**: peak-to-peak values below 10 uV or above 1000 uV discarded
- **Excluded**: A/B 100 ms (only 1 subject)
- **41a73 block correction**: R config updated to match MATLAB block definitions (2 conditioning protocols instead of 4)

### Model

Linear mixed-effects model on log-transformed cell-median EP amplitude (log uV):

```
PPvec ~ mapStimLevel + chanInCond + overallBlockType * pre_post
        + (1 + stim_linpoly|subjectNum:chanVec) + (1|subjectNum) + (1|subjectNum:blockVec)
```

Crossed random-effect structure: channels and blocks are crossed within subjects because multi-channel subjects had both channels recorded simultaneously during each block. Random linear dose-response slope per channel allows each electrode to have its own dose-response steepness, reflecting cortical position. `stim_linpoly` is the linear polynomial contrast extracted from `contr.poly()`. Model comparison: channel slope (AIC=207) beats block slope (AIC=277) and intercept-only (AIC=347). Log transformation is required for the crossed RE to converge (raw-scale model is singular: block RE variance = 0).

### Random Effects

| Component | Variance | SD | Corr | Groups |
|-----------|----------|-----|------|--------|
| Subject:Block (intercept) | 0.047 | 0.217 | | 51 |
| Subject:Channel (intercept) | 0.405 | 0.636 | | 17 |
| Subject:Channel (slope) | 0.323 | 0.569 | -0.35 | |
| Subject | 0.117 | 0.342 | | 13 |
| Residual | 0.049 | 0.221 | | -- |
| **Total** | **0.941** | **0.970** | | -- |

Model is not singular. Shapiro-Wilk on residuals: W=0.991, p=0.091 (passes normality). Skewness=-0.18, excess kurtosis=0.34.

### Type III ANOVA (Satterthwaite df)

| Effect | Sum Sq | Mean Sq | NumDF | DenDF | F | p |
|--------|--------|---------|-------|-------|---|---|
| Stimulation level | 3.175 | 1.058 | 3 | 34.3 | 21.8 | 4.5e-08 |
| Channel in conditioning pair | 0.362 | 0.362 | 1 | 9.8 | 7.45 | 0.022 |
| Conditioning protocol | 0.441 | 0.147 | 3 | 41.0 | 3.02 | 0.041 |
| Pre/post | 0.060 | 0.060 | 1 | 34.4 | 1.22 | 0.276 |
| Protocol x pre/post | 0.208 | 0.069 | 3 | 37.8 | 1.42 | 0.251 |

### Fixed Effects

Sum-to-zero coding (`contr.sum`): intercept is the grand mean; each coefficient is the deviation of that level from the grand mean. The last level of each factor is implicit (negative sum of the others).

| Predictor | Estimate | SE | df | t | p |
|-----------|----------|-----|-----|---|---|
| Intercept (grand mean) | 4.681 | 0.200 | 14.0 | 23.44 | < 1e-10 |
| Stimulation level (linear) | 0.954 | 0.143 | 13.9 | 6.66 | 1.1 x 10^-05 |
| Stimulation level (quadratic) | -0.128 | 0.028 | 190.6 | -4.61 | 7.2 x 10^-06 |
| Stimulation level (cubic) | -0.0001 | 0.027 | 189.2 | -0.004 | 0.997 |
| Channel not in pair | -1.135 | 0.416 | 9.8 | -2.73 | 0.022 |
| A/A 200 (dev. from grand mean) | 0.142 | 0.085 | 36.4 | 1.68 | 0.101 |
| A/A 25 (dev. from grand mean) | -0.101 | 0.080 | 45.6 | -1.26 | 0.215 |
| A/B 200 (dev. from grand mean) | 0.137 | 0.061 | 35.5 | 2.24 | 0.032 |
| Post (dev. from grand mean) | 0.041 | 0.037 | 34.4 | 1.11 | 0.276 |
| A/A 200 x post | -0.046 | 0.065 | 34.1 | -0.71 | 0.483 |
| A/A 25 x post | 0.059 | 0.072 | 43.8 | 0.82 | 0.418 |
| A/B 200 x post | 0.066 | 0.053 | 33.6 | 1.24 | 0.225 |

A/B 25 coefficients are implicit: condition deviation = -(0.142 + (-0.101) + 0.137) = -0.178; interaction = -((-0.046) + 0.059 + 0.066) = -0.079.

### Pre-Post Contrasts Within Each Condition (Satterthwaite df)

| Condition | Estimate (post - pre) | SE | 95% CI | df | Cohen's d | d 95% CI |
|-----------|----------------------|-----|---------|-----|-----------|----------|
| A/A 200 | -0.011 | 0.152 | [-0.320, 0.299] | 33.9 | -0.013 | [-0.392, 0.367] |
| A/A 25 | 0.200 | 0.183 | [-0.169, 0.569] | 41.9 | 0.239 | [-0.218, 0.696] |
| A/B 200 | 0.213 | 0.107 | [-0.004, 0.431] | 32.8 | 0.255 | [-0.024, 0.534] |
| A/B 25 | -0.074 | 0.119 | [-0.316, 0.168] | 37.5 | -0.088 | [-0.389, 0.212] |

All estimates and effect sizes are on the log(uV) scale. Cohen's d uses marginal total SD (0.835) as denominator, with proper E[x^2] weighting for the random slope variance component. All CIs cross zero.

### Conditioning Protocol Pairwise Comparisons (Tukey-adjusted)

| Comparison | Estimate | SE | df | t | p |
|------------|----------|-----|-----|---|---|
| A/A 200 - A/A 25 | 0.243 | 0.146 | 39.5 | 1.66 | 0.358 |
| A/A 200 - A/B 200 | 0.006 | 0.104 | 34.8 | 0.05 | 1.000 |
| A/A 200 - A/B 25 | 0.320 | 0.136 | 38.7 | 2.35 | 0.105 |
| A/A 25 - A/B 200 | -0.237 | 0.122 | 40.8 | -1.95 | 0.224 |
| A/A 25 - A/B 25 | 0.077 | 0.108 | 55.4 | 0.71 | 0.893 |
| A/B 200 - A/B 25 | 0.314 | 0.107 | 39.0 | 2.93 | 0.028 |

### Degrees of Freedom by Hierarchy Level

| Level | Effect | Satterthwaite df | Rationale |
|-------|--------|-----------------|-----------|
| Within-block | Stimulation level (Q, C) | ~190 | 262 cell medians - 51 block groups - params |
| Channel (stim slope) | Stimulation level (L) | ~14 | 17 channels with random slope for linear dose |
| Block (within-subject) | Conditioning protocol, pre/post, interaction | ~34-41 | 51 unique sessions minus params |
| Channel (within-subject) | Channel in conditioning pair | ~10 | 17 channels |

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
| 1 | -28.9 | 0.166 | 0.429 |
| 2 | 27.0 | 0.429 | 0.429 |
| 3 | -28.1 | 0.322 | 0.429 |
| 4 | 11.9 | 0.301 | 0.429 |

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
| 5-min conditioning (stratified pooled) | -4.5 | 0.719 |
| 15-min conditioning (stratified pooled) | 64.4 | < 0.0001 |
| Interaction: 15-min minus 5-min (stratified pooled) | 68.9 | 0.0001 |

### Interpretation

15-minute conditioning produces a large, robust increase in EP amplitude (~82-86 uV at stim levels 2-4, p < 0.007; all three survive BH-FDR within the family of 4 stim levels, p_BH <= 0.008). 5-minute conditioning produces **no significant effect** (stratified pooled = -4.5 uV, p = 0.719; per-stim-level effects are mixed in direction and all non-significant before and after BH-FDR). The stratified pooled interaction confirms that 15-min conditioning produces a significantly larger effect than 5-min (68.9 uV, p = 0.0001). At the per-stim-level, the interaction is robust at stim level 3 (p_BH = 0.004) and marginal at stim level 4 (raw p = 0.042, p_BH = 0.084); the stratified pooled test is the primary inference.

**Note on data provenance:** The previous `_avg_5.csv` data files had unknown provenance (no traceable MATLAB code generated them). When regenerated from the current committed code (with proper bad trial exclusion and known Savitzky-Golay parameters), the 5-minute conditioning effect — previously reported as significant (p=0.004) — became non-significant (p=0.719). The 15-minute effect and interaction remained robust. This underscores the importance of reproducible data extraction pipelines.

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
