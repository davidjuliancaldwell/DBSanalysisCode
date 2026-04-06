# Statistical Results Summary: DBS Paired Pulse EP Analysis

## 1. Primary Multi-Subject Analysis

### Study Design

- **Subjects**: N=9 (after filtering; 12 total enrolled, 8 PD, 1 MD in final sample)
- **Channels**: 13 subject-channel combinations (9 subjects with 1 channel, 4 with 2 channels recorded simultaneously)
- **Recording sessions**: 37 unique subject-block sessions
- **Conditions**: 4 paired-pulse conditioning protocols
  - A/A 200 ms (same-polarity control, long ISI): 2 subjects, 4 blocks
  - A/A 25 ms (same-polarity control, short ISI): 3 subjects, 6 blocks
  - A/B 200 ms (opposite-polarity, long ISI): 8 subjects, 16 blocks
  - A/B 25 ms (opposite-polarity, short ISI): 5 subjects, 12 blocks
- **Observations**: 124 cell medians aggregated from 1,462 trial-level observations (5-trial sequential averages, no trimming)
- **Stimulation levels**: 2 ordered levels per block (after min_stim_level=3 filter)
- **Excluded**: A/B 100 ms (only 1 subject)

### Model

Linear mixed-effects model on log-transformed cell-median EP amplitude (log uV):

```
PPvec ~ mapStimLevel + chanInCond + overallBlockType * pre_post
        + (1|subjectNum/chanVec) + (1|subjectNum:blockVec)
```

Crossed random-effect structure: channels and blocks are crossed within subjects because multi-channel subjects had both channels recorded simultaneously during each block. Log transformation is required for the crossed RE to converge (raw-scale model is singular: block RE variance = 0).

### Random Effects

| Component | Variance | SD | Groups |
|-----------|----------|-----|--------|
| Subject | 0.183 | 0.428 | 9 |
| Subject:Channel | 0.377 | 0.614 | 13 |
| Subject:Block | 0.031 | 0.176 | 45 |
| Residual | 0.055 | 0.236 | -- |
| **Total** | **0.646** | **0.804** | -- |

Model is not singular. ICC (subject + channel + block) = 0.914.

### Type III ANOVA (Satterthwaite df)

| Effect | Sum Sq | F | NumDF | DenDF | p |
|--------|--------|---|-------|-------|---|
| Stimulation level | 2.839 | 51.18 | 1 | 77.8 | 4.0 x 10^-10 |
| Channel in conditioning pair | 0.462 | 8.33 | 1 | 6.9 | 0.024 |
| Conditioning protocol | 0.406 | 2.44 | 3 | 38.8 | 0.079 |
| Pre/post | 0.032 | 0.57 | 1 | 33.2 | 0.455 |
| Protocol x pre/post | 0.431 | 2.59 | 3 | 36.1 | 0.068 |

### Fixed Effects

| Predictor | Estimate | SE | df | t | p |
|-----------|----------|-----|-----|---|---|
| Intercept | 5.203 | 0.274 | 14.9 | 18.99 | 7.6 x 10^-12 |
| Stimulation level (linear) | 0.214 | 0.030 | 77.8 | 7.15 | 4.1 x 10^-10 |
| Channel not in pair | -1.284 | 0.445 | 6.9 | -2.89 | 0.024 |
| A/A 25 (vs A/A 200) | -0.078 | 0.192 | 34.8 | -0.40 | 0.689 |
| A/B 200 (vs A/A 200) | 0.116 | 0.147 | 34.6 | 0.79 | 0.436 |
| A/B 25 (vs A/A 200) | -0.329 | 0.160 | 36.1 | -2.06 | 0.047 |
| Pre (vs post) | -0.069 | 0.165 | 34.7 | -0.42 | 0.680 |
| A/A 25 x pre | -0.155 | 0.247 | 36.9 | -0.63 | 0.535 |
| A/B 200 x pre | -0.101 | 0.196 | 33.6 | -0.52 | 0.610 |
| A/B 25 x pre | 0.303 | 0.206 | 35.1 | 1.47 | 0.150 |

### Pre-Post Contrasts Within Each Condition (Kenward-Roger df)

| Condition | Estimate (post - pre) | SE | 95% CI | df | t | p | Cohen's d | d 95% CI |
|-----------|----------------------|-----|---------|-----|---|---|-----------|----------|
| A/A 200 | 0.069 | 0.165 | [-0.268, 0.406] | 30.2 | 0.42 | 0.681 | 0.085 | [-0.359, 0.529] |
| A/A 25 | 0.223 | 0.186 | [-0.155, 0.601] | 34.1 | 1.20 | 0.238 | 0.278 | [-0.223, 0.779] |
| A/B 200 | 0.170 | 0.107 | [-0.050, 0.390] | 26.9 | 1.59 | 0.125 | 0.211 | [-0.090, 0.513] |
| A/B 25 | -0.234 | 0.124 | [-0.486, 0.018] | 31.4 | -1.89 | 0.068 | -0.291 | [-0.637, 0.054] |

All estimates and effect sizes are on the log(uV) scale. Cohen's d uses total SD (0.804) as denominator for cross-study comparability. All CIs cross zero.

### Conditioning Protocol Pairwise Comparisons (Tukey-adjusted)

| Comparison | Estimate | SE | df | t | p |
|------------|----------|-----|-----|---|---|
| A/A 200 - A/A 25 | 0.155 | 0.145 | 35.0 | 1.07 | 0.710 |
| A/A 200 - A/B 200 | -0.065 | 0.110 | 30.8 | -0.59 | 0.933 |
| A/A 200 - A/B 25 | 0.178 | 0.123 | 33.3 | 1.45 | 0.477 |
| A/A 25 - A/B 200 | -0.220 | 0.121 | 34.6 | -1.83 | 0.279 |
| A/A 25 - A/B 25 | 0.023 | 0.114 | 40.8 | 0.20 | 0.997 |
| A/B 200 - A/B 25 | 0.243 | 0.094 | 31.9 | 2.58 | 0.067 |

### Degrees of Freedom by Hierarchy Level

| Level | Effect | Satterthwaite df | Rationale |
|-------|--------|-----------------|-----------|
| Within-block | Stimulation level | ~78 | 124 cell medians - 45 block groups - params |
| Block (within-subject) | Conditioning protocol, pre/post, interaction | ~33-39 | 37 unique sessions minus params |
| Channel (within-subject) | Channel in conditioning pair | ~7 | 13 channels, only 3 subjects with both levels |

---

## 2. Single-Subject Conditioning Length (subject a23ed)

### Study Design

- **Subject**: 1 (a23ed)
- **Channel**: 1 (channel 5)
- **Conditions**: 5-min vs 15-min A/B 200 ms conditioning, each with its own baseline
  - Block 2 (baseline) vs block 3 (post 5-min conditioning)
  - Block 5 (baseline) vs block 6 (post 15-min conditioning)
- **Trials per cell**: 8-10 (after 5-trial averaging, no trimming)
- **Stimulation levels**: 3 (levels 2, 3, 4)
- **Total trials**: 110
- **Data scale**: raw uV (not log-transformed)

### Permutation Tests (Primary Analysis)

10,000 permutations per test. Test statistic: difference of medians. Block labels shuffled within each stimulation level. Two-sided p-values: `mean(abs(perm) >= abs(obs))`. Trial-level permutation is the only feasible approach (block-level gives C(2,1) = 2 arrangements per comparison, minimum p = 1.0).

**Per-stim-level tests (permutation within each stim level):**

**5-minute conditioning (post - baseline):**

| Stim level | Observed diff (uV) | p |
|------------|-------------------|---|
| 2 | 22.0 | 0.055 |
| 3 | 8.3 | 0.749 |
| 4 | 19.3 | 0.045 |

**15-minute conditioning (post - baseline):**

| Stim level | Observed diff (uV) | p |
|------------|-------------------|---|
| 2 | 78.8 | < 0.001 |
| 3 | 84.2 | < 0.001 |
| 4 | 68.5 | 0.001 |

**Per-stim-level interaction: (15-min effect) - (5-min effect):**

| Stim level | Observed interaction (uV) | p |
|------------|--------------------------|---|
| 2 | 56.8 | 0.053 |
| 3 | 75.9 | 0.135 |
| 4 | 49.2 | 0.176 |

**Stratified pooled tests:** These pool across stim levels while preserving stim-level balance. Within each permutation iteration, labels are shuffled independently within each stim level, per-stim-level median differences are computed, then averaged. This avoids the problem with naive pooling, where random stim-level imbalance widens the null distribution and makes the test conservative.

| Test | Observed (uV) | p |
|------|---------------|------|
| 5-min conditioning (stratified pooled) | 16.5 | 0.013 |
| 15-min conditioning (stratified pooled) | 77.2 | < 0.0001 |
| Interaction: 15-min minus 5-min (stratified pooled) | 60.7 | 0.0034 |

### Interpretation

15-minute conditioning produces a large, robust increase in EP amplitude (~70-85 uV, p < 0.001 at all stimulation levels). 5-minute conditioning produces a smaller effect that is significant only at the highest stimulation level per-stim-level (p=0.045) but reaches significance in the stratified pooled test (p=0.013). The stratified pooled interaction confirms that 15-min conditioning produces a significantly larger effect than 5-min (p=0.003), whereas the per-stim-level interaction tests are underpowered (~10 trials per cell).

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
- **Total trials**: 356 (after 5-trial averaging, filtering stim level >= 2)
- **Stim levels**: 3 (levels 2, 3, 4)
- **Trials per (channel, stim level)**: ~58-60 (22-24 awake, 35-36 asleep)
- **Data scale**: raw uV (`log_data = FALSE`; log produces skewness=-5.7 and kurtosis=65 for this subject due to near-zero EP values)
- **No trimming** (`trim_data = FALSE`)

### LMM (Model-Based Analysis)

```
fit.lmm1 = PPvec ~ mapStimLevel + blockType + chanVec + (1|blockVec)  # additive (primary)
fit.lmm2 = PPvec ~ mapStimLevel + blockType * mapStimLevel + chanVec + (1|blockVec)  # interaction
```

Trial-level data with block RE (10 groups). `chanVec` as fixed effect (only 2 channels). Blocks are perfectly nested within channels (odd blocks = chan 6, even = chan 4); unique block IDs make `(1|blockVec)` implicitly handle this. Additive model preferred over interaction (LRT p = 0.89, AIC/BIC favor additive). Not singular. Random effects: block SD = 13.2, residual SD = 46.5, ICC = 0.074.

**ANOVA (Type III, Satterthwaite):**

| Effect | F | NumDF | DenDF | p |
|--------|---|-------|-------|---|
| mapStimLevel | 290.4 | 2 | 344 | < 2e-16 |
| blockType (awake vs asleep) | 0.48 | 1 | 7.0 | 0.51 |
| chanVec | 64.7 | 1 | 7.0 | 8.9e-05 |

**Degrees of freedom:** Satterthwaite correctly assigns ~7 df for between-block effects (blockType, chanVec) and ~344 df for within-block effects (mapStimLevel). The between-block df = 10 blocks - 3 between-block params (intercept + blockType + chanVec) = 7. ICC > 0 confirms the RE is estimated successfully and df are not artificially inflated.

Awake vs asleep (blockType): estimate = 6.8 uV, SE = 9.9, df = 7.0, t = 0.69, p = 0.51. **Not significant.**

### Permutation Tests (Trial-Level, Stratified by Stim Level)

Trial-level permutation tests per channel (channels are independent). 10,000 permutations. Test statistic: difference of medians (awake - asleep). Awake/asleep labels shuffled within each (channel, stim_level) stratum. Two-sided p-values: `mean(abs(perm) >= abs(obs))`.

**Per channel, per stim level:**

| Channel | Stim Level | Observed (uV) | p |
|---------|-----------|---------------|------|
| 6 | 2 | 12.7 | 0.52 |
| 6 | 3 | 14.3 | 0.38 |
| 6 | 4 | 18.6 | 0.27 |
| 4 | 2 | -0.8 | 0.94 |
| 4 | 3 | 12.7 | 0.25 |
| 4 | 4 | 5.9 | 0.31 |

**Stratified pooled (per channel):** Independently permutes within each stim level, averages per-stim-level median differences. Preserves stim-level balance in every permutation.

| Channel | Observed (uV) | p |
|---------|---------------|------|
| 6 | 15.2 | 0.12 |
| 4 | 5.9 | 0.21 |

**Consistent with LMM:** No significant awake/asleep effect. Channel 6 shows a non-significant trend. Trial-level permutation is mildly liberal due to ICC ~ 0.07 (effective n inflation ~1.77x), but since results are clearly non-significant, this does not affect interpretation. Block-level permutation has insufficient resolution (C(5,2) = 10 per channel, minimum p = 0.10).

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
| 15-min conditioning increases EPs in subject a23ed | < 0.001 |
| No individual conditioning protocol reaches significance (multi-subject) | > 0.05 |

### Findings Sensitive to Analytic Choices

| Finding | Range of p across configurations | Status |
|---------|--------------------------------|--------|
| Conditioning protocol x pre/post interaction (multi-subject) | 0.05 - 0.45 | Exploratory |
| 15-min vs 5-min interaction (a23ed) | 0.003 - 0.18 | Stratified pooled p=0.003; per-stim-level p=0.05-0.18. Depends on pooling strategy and stim level |
| No awake/asleep effect (3d413) | 0.12 - 0.94 | Consistent across LMM (p=0.51) and permutation (per-channel pooled p=0.12, 0.21) |

### Key Analytic Decisions and Their Impact

| Decision | Alternatives tested | Impact |
|----------|-------------------|--------|
| Cell median vs cell mean | Both | Minimal for main effects; marginal interaction shifts ~0.01 in p |
| Trim vs no trim | Both | Adds/removes ~20% of cell medians due to n>=10 filter; minimal impact on conclusions |
| Log vs raw scale (main) | Both | Raw model is singular (block RE = 0); log required for proper RE convergence |
| Log vs raw scale (3d413) | Both | Log produces catastrophic non-normality (kurtosis=65); raw is correct for single-subject |
| 3 vs 4 conditions | Both | 4 conditions (with A/A 25) uses 3 df for interaction vs 2; dilutes F-statistic |
| Flat vs crossed RE | Both | Flat gives anti-conservative df for channel effects (37 vs 7); crossed is correct |
| Permutation: mean vs median | Both | 15-min effect robust to both; interaction significant only with mean |
| Permutation: stratified vs naive pooling | Both | Naive pooling mixes stim levels, widening null (conservative). Stratified pooling preserves stim-level balance, giving tighter null and cleaner inference. Stratified pooled interaction p=0.003 vs naive approach |
| Trial-level vs block-level permutation (3d413) | Both | Block-level: C(5,2)=10 per channel, minimum p=0.10 (too coarse). Trial-level: adequate resolution but mildly liberal (ICC~0.07). Trial-level chosen |
