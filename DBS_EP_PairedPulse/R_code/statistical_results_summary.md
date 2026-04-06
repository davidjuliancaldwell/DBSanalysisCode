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

10,000 permutations per test. Test statistic: difference of medians. Block labels shuffled within each stimulation level. Two-sided p-values.

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

**Interaction: (15-min effect) - (5-min effect):**

| Stim level | Observed interaction (uV) | p |
|------------|--------------------------|---|
| 2 | 56.8 | 0.053 |
| 3 | 75.9 | 0.135 |
| 4 | 49.2 | 0.176 |

### Interpretation

15-minute conditioning produces a large, robust increase in EP amplitude (~70-85 uV, p < 0.001 at all stimulation levels). 5-minute conditioning produces a smaller effect that is significant only at the highest stimulation level (p=0.045). The per-stim-level interaction tests (whether 15-min is larger than 5-min) are marginal to non-significant, limited by the small number of trials per cell (~10) available for the median statistic.

### Caveats

- Permutation tests shuffle individual trials between blocks, assuming trial exchangeability. Within-block temporal correlation (partially mitigated by 5-trial averaging) could make p-values anti-conservative.
- With 4 blocks, mixed-effects models have convergence issues. The permutation approach avoids model-based assumptions but shares the pseudoreplication concern.
- Single-subject results do not generalize without replication.

---

## 3. Single-Subject Anesthesia Variability (subject 3d413)

### Study Design

- **Subject**: 1 (3d413)
- **Channels**: 2 (channels 4 and 6)
- **Blocks**: 10 (5 per channel: 3 asleep, 2 awake)
- **Data scale**: raw uV (`log_data = FALSE`; log produces skewness=-5.7 and kurtosis=65 for this subject due to near-zero EP values)
- **No trimming** (`trim_data = FALSE`)

### Model

```
PPvec ~ mapStimLevel + blockType + chanVec + (1|blockVec)
```

Trial-level data with block RE (10 groups -- sufficient for RE estimation). `chanVec` as fixed effect (only 2 channels). Additive model preferred over interaction model (blockType x mapStimLevel) by AIC and BIC. Not singular. Residuals: skewness=-0.5, excess kurtosis=2.9.

Awake vs asleep (blockType): estimate=6.8 uV, SE=9.9, df=7.0, p=0.51. Not significant.

### Effect Sizes

Cohen's d (awake vs asleep) computed per (channel, stim level) on trial-level data. Raw uV scale.

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
| 15-min vs 5-min interaction (a23ed) | 0.003 - 0.18 | Depends on statistic (mean vs median) and per-stim-level vs pooled |

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
