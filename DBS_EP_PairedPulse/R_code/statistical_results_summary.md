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
- **Observations**: 250 cell medians (5-trial sequential averages, no trimming, n>=3 per cell filter, regenerated data)
- **Stimulation levels**: all available ordered levels included (min_stim_level=1, up to 4 per subject)
- **Amplitude filter**: peak-to-peak values below 10 uV or above 1000 uV discarded
- **Excluded**: A/B 100 ms (only 1 subject)
- **41a73 block correction**: R config updated to match MATLAB block definitions (2 conditioning protocols instead of 4)

### Model

Linear mixed-effects model on log-transformed cell-median EP amplitude (log uV):

```
PPvec ~ mapStimLevel + chanInCond + overallBlockType * pre_post
        + (1 + stim_linpoly|subjectNum:chanVec) + (1|subjectNum) + (1|subjectNum:blockVec)
```

Crossed random-effect structure: channels and blocks are crossed within subjects because multi-channel subjects had both channels recorded simultaneously during each block. Random linear dose-response slope per channel allows each electrode to have its own dose-response steepness, reflecting cortical position. `stim_linpoly` is the linear polynomial contrast extracted from `contr.poly()`. Model comparison: channel slope (AIC=158) beats block slope (AIC=210) and intercept-only (AIC=246). When both slopes are included (AIC=161), the block slope collapses to SD~0.04, confirming dose-response steepness is a channel property, not a session property. Log transformation is required for the crossed RE to converge (raw-scale model is singular: block RE variance = 0).

### Random Effects

| Component | Variance | SD | Corr | Groups |
|-----------|----------|-----|------|--------|
| Subject:Block (intercept) | 0.047 | 0.217 | | 51 |
| Subject:Channel (intercept) | 0.359 | 0.599 | | 17 |
| Subject:Channel (slope) | 0.602 | 0.776 | -0.11 | |
| Subject | 0.111 | 0.333 | | 13 |
| Residual | 0.055 | 0.234 | | -- |
| **Total** | **1.174** | **1.083** | | -- |

Model is not singular. Shapiro-Wilk on residuals: W=0.993, p=0.259 (passes normality). Skewness=-0.22, excess kurtosis=0.24.

### Type III ANOVA (Satterthwaite df)

| Effect | Sum Sq | Mean Sq | NumDF | DenDF | F | p |
|--------|--------|---------|-------|-------|---|---|
| Stimulation level | 2.931 | 0.977 | 3 | 35.5 | 17.8 | 3.2e-07 |
| Channel in conditioning pair | 0.453 | 0.453 | 1 | 9.9 | 8.25 | 0.017 |
| Conditioning protocol | 0.485 | 0.162 | 3 | 41.8 | 2.95 | 0.044 |
| Pre/post | 0.069 | 0.069 | 1 | 35.0 | 1.27 | 0.268 |
| Protocol x pre/post | 0.234 | 0.078 | 3 | 38.5 | 1.42 | 0.251 |

### Fixed Effects

| Predictor | Estimate | SE | df | t | p |
|-----------|----------|-----|-----|---|---|
| Intercept | 4.736 | 0.232 | 27.3 | 20.41 | < 1e-10 |
| Stimulation level (linear) | 0.878 | 0.143 | 14.6 | 6.15 | 2.4 x 10^-05 |
| Stimulation level (quadratic) | -0.067 | 0.031 | 187.9 | -2.17 | 0.031 |
| Stimulation level (cubic) | 0.081 | 0.031 | 186.7 | 2.56 | 0.011 |
| Channel not in pair | -1.130 | 0.393 | 9.9 | -2.87 | 0.017 |
| A/A 25 (vs A/A 200) | -0.194 | 0.192 | 36.4 | -1.01 | 0.320 |
| A/B 200 (vs A/A 200) | 0.072 | 0.139 | 37.1 | 0.52 | 0.609 |
| A/B 25 (vs A/A 200) | -0.398 | 0.170 | 38.2 | -2.35 | 0.024 |
| Pre (vs post) | 0.003 | 0.156 | 35.6 | 0.02 | 0.987 |
| A/A 25 x pre | -0.206 | 0.242 | 41.0 | -0.85 | 0.401 |
| A/B 200 x pre | -0.219 | 0.189 | 35.5 | -1.16 | 0.254 |
| A/B 25 x pre | 0.076 | 0.196 | 37.4 | 0.39 | 0.701 |

### Pre-Post Contrasts Within Each Condition (Satterthwaite df)

| Condition | Estimate (post - pre) | SE | 95% CI | df | Cohen's d | d 95% CI |
|-----------|----------------------|-----|---------|-----|-----------|----------|
| A/A 200 | -0.003 | 0.156 | [-0.319, 0.314] | 35.6 | -0.002 | [-0.301, 0.297] |
| A/A 25 | 0.203 | 0.186 | [-0.172, 0.577] | 41.8 | 0.19 | [-0.170, 0.545] |
| A/B 200 | 0.217 | 0.109 | [-0.004, 0.438] | 33.6 | 0.20 | [-0.018, 0.418] |
| A/B 25 | -0.078 | 0.121 | [-0.323, 0.167] | 37.7 | -0.07 | [-0.307, 0.163] |

All estimates and effect sizes are on the log(uV) scale. Cohen's d uses total SD (1.083) as denominator for cross-study comparability. All CIs cross zero.

### Conditioning Protocol Pairwise Comparisons (Tukey-adjusted)

| Comparison | Estimate | SE | df | t | p |
|------------|----------|-----|-----|---|---|
| A/A 200 - A/A 25 | 0.241 | 0.149 | 40.3 | 1.62 | 0.381 |
| A/A 200 - A/B 200 | 0.002 | 0.106 | 36.5 | 0.02 | 1.000 |
| A/A 200 - A/B 25 | 0.316 | 0.138 | 39.8 | 2.29 | 0.119 |
| A/A 25 - A/B 200 | -0.238 | 0.124 | 41.0 | -1.93 | 0.232 |
| A/A 25 - A/B 25 | 0.076 | 0.111 | 55.2 | 0.69 | 0.902 |
| A/B 200 - A/B 25 | 0.314 | 0.109 | 39.3 | 2.90 | 0.030 |

### Degrees of Freedom by Hierarchy Level

| Level | Effect | Satterthwaite df | Rationale |
|-------|--------|-----------------|-----------|
| Within-block | Stimulation level (Q, C) | ~187 | 250 cell medians - 51 block groups - params |
| Channel (stim slope) | Stimulation level (L) | ~15 | 17 channels with random slope for linear dose |
| Block (within-subject) | Conditioning protocol, pre/post, interaction | ~35-42 | 51 unique sessions minus params |
| Channel (within-subject) | Channel in conditioning pair | ~10 | 17 channels |

### Secondary Model: halfBlock (Temporal Modulation Within Blocks)

Tests whether the conditioning effect differs between the first and second half of trials within each block. Trials within each (block, stim level) cell are split at the midpoint; cell medians computed per half. Uses `dataListAggHalf` (500 obs from 250 cells × 2 halves).

```
PPvec ~ mapStimLevel + chanInCond + overallBlockType * pre_post * halfBlock
        + (1 + stim_linpoly|subjectNum:chanVec) + (1|subjectNum) + (1|subjectNum:blockVec)
```

| Effect | F | NumDF | DenDF | p |
|--------|---|-------|-------|---|
| halfBlock | 0.01 | 1 | 417 | 0.909 |
| overallBlockType:halfBlock | 0.22 | 3 | 417 | 0.879 |
| pre_post:halfBlock | 0.67 | 1 | 417 | 0.414 |
| overallBlockType:pre_post:halfBlock | 1.14 | 3 | 417 | 0.334 |

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

| Stim level | Observed diff (uV) | p |
|------------|-------------------|---|
| 1 | -28.9 | 0.166 |
| 2 | 27.0 | 0.429 |
| 3 | -28.1 | 0.322 |
| 4 | 11.9 | 0.301 |

**15-minute conditioning (post - baseline):**

| Stim level | Observed diff (uV) | p |
|------------|-------------------|---|
| 1 | 4.4 | 0.541 |
| 2 | 84.4 | < 0.001 |
| 3 | 82.3 | < 0.001 |
| 4 | 86.4 | 0.006 |

**Per-stim-level interaction: (15-min effect) - (5-min effect):**

| Stim level | Observed interaction (uV) | p |
|------------|--------------------------|---|
| 1 | 33.3 | 0.124 |
| 2 | 57.5 | 0.228 |
| 3 | 110.4 | 0.001 |
| 4 | 74.5 | 0.042 |

**Stratified pooled tests:** These pool across stim levels while preserving stim-level balance. Within each permutation iteration, labels are shuffled independently within each stim level, per-stim-level median differences are computed, then averaged. This avoids the problem with naive pooling, where random stim-level imbalance widens the null distribution and makes the test conservative.

| Test | Observed (uV) | p |
|------|---------------|------|
| 5-min conditioning (stratified pooled) | -4.5 | 0.719 |
| 15-min conditioning (stratified pooled) | 64.4 | < 0.0001 |
| Interaction: 15-min minus 5-min (stratified pooled) | 68.9 | 0.0001 |

### Interpretation

15-minute conditioning produces a large, robust increase in EP amplitude (~82-86 uV at stim levels 2-4, p < 0.007). 5-minute conditioning produces **no significant effect** (stratified pooled = -4.5 uV, p = 0.719; per-stim-level effects are mixed in direction and all non-significant). The stratified pooled interaction confirms that 15-min conditioning produces a significantly larger effect than 5-min (68.9 uV, p = 0.0001).

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

Awake vs asleep (blockType): estimate = 6.8 uV, SE = 8.2, df = 6.6, t = 0.83, p = 0.436. **Not significant.**

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

All scripts use R defaults: `as.ordered()` → polynomial contrasts (contr.poly), `as.factor()` → treatment contrasts (contr.treatment). No scripts override global contrast settings. Contrasts are per-variable in R — no resetting needed between variables.

- **`mapStimLevel`** (`as.ordered`): polynomial contrasts (.L linear, .Q quadratic). Appropriate for ordinal dose-response. Type III F-test is contrast-invariant; only individual coefficient estimates (e.g., the linear trend) use polynomial decomposition.
- **All other fixed-effect factors** (`as.factor`): treatment (dummy) contrasts. All 2-level factors (blockType, chanVec, pre_post, chanInCond) have 1 df — contrast type is irrelevant for the F-test.
- **RE grouping factors** (subjectNum, blockVec): contrasts are unused — they only define grouping structure for random intercepts.

---

## Robustness Assessment

### Findings Robust Across All Analytic Configurations

| Finding | Minimum p across configurations |
|---------|-------------------------------|
| Stimulation level drives EP amplitude | < 10^-6 |
| Channels in conditioning pair have larger EPs | ~ 0.02 |
| A/B 200 vs A/B 25 overall amplitude differs | 0.030 (Tukey) |
| 15-min conditioning increases EPs in subject a23ed | < 0.001 (stim levels 2-4) |
| No individual conditioning protocol pre-post change reaches significance | > 0.05 (A/B 200 approaches: d=0.20, CI [-0.018, 0.418]) |
| No significant conditioning x pre/post interaction | 0.251 |
| No temporal modulation within blocks (halfBlock) | Three-way interaction p=0.334; halfBlock main p=0.909 |

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
| All stim levels vs min_stim=3 | Both | Including all levels increases obs (124→202) and adds cubic polynomial term; conclusions unchanged |
| 41a73 block correction | Before/after | R config corrected to match MATLAB; removes phantom A/B 25 and A/A 25 conditions for this subject |
| Log vs raw scale (main) | Both | Raw model is singular (block RE = 0); log required for proper RE convergence |
| Log vs raw scale (3d413) | Both | Log produces catastrophic non-normality (kurtosis=65); raw is correct for single-subject |
| 3 vs 4 conditions | Both | 4 conditions (with A/A 25) uses 3 df for interaction vs 2; dilutes F-statistic |
| Flat vs crossed RE | Both | Flat gives anti-conservative df for channel effects (37 vs 7); crossed is correct |
| Permutation: mean vs median | Both | 15-min effect robust to both; interaction significant only with mean |
| Permutation: stratified vs naive pooling | Both | Naive pooling mixes stim levels, widening null (conservative). Stratified pooling preserves stim-level balance, giving tighter null and cleaner inference. Stratified pooled interaction p=0.003 vs naive approach |
| Trial-level vs block-level permutation (3d413) | Both | Block-level: C(5,2)=10 per channel, minimum p=0.10 (too coarse). Trial-level: adequate resolution but mildly liberal (ICC~0.07). Trial-level chosen |
