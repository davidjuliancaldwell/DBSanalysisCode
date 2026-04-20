# Permutation tests (median post - pre) for every (subject, channel, pre/post-pair, stim-level)
# combination described by the per-subject R config files. BH-FDR applied within each
# (subject, channel, condition) family of 4 stim levels.
#
# Mirrors the per-stim-level test used in length_conditioning_a23ed.R:
#   obs_diff  = median(post_trials) - median(pre_trials)
#   perm_diff = same, with block labels shuffled within the stim-level stratum
#   two-sided p = mean(|perm| >= |obs|)
#
# Trial-level data = 5-trial sequential averages (PPvec column of
# `_PairedPulseData_avg_5.csv`). Amplitude filter 10-1000 uV applied first.
# Subjects / channels / compare pairs come from R_config_files/subj_{sid}.R.
# 3d413 is excluded (awake/asleep study, no conditioning) to match the
# main analysis sidVec.
#
# Output:
#   R_output/permutation_results_all_subjects.csv
#   R_output/permutation_results_all_subjects_stim4.csv (subset, stim level 4)
#
# Run time: ~3-5 minutes with nPerm = 10000.

library(here)
library(dplyr)

set.seed(42)
nPerm <- 10000

# Matches dose_response_R_script_trim_conditions.R sidVec (excludes 3d413).
sidVec <- c('46c2a','c963f','2e114','fe7df','e6f3c','9f852',
            '8e907','08b13','e9c9b','41a73','68574',
            '01fee','a23ed')

dataDir <- here("DBS_EP_PairedPulse","R_data")
configDir <- here("DBS_EP_PairedPulse","R_config_files")
outDir <- here("DBS_EP_PairedPulse","R_output")
if (!dir.exists(outDir)) dir.create(outDir, recursive = TRUE)

all_results <- list()

for (sid in sidVec) {

  configFile <- file.path(configDir, paste0("subj_", sid, ".R"))
  if (!file.exists(configFile)) {
    cat(sprintf("[skip] %s: no subj config\n", sid))
    next
  }
  rm(list = intersect(c("subjectNum","chanIntVec","blockIntLM","blockType","whichCompareVec"),
                      ls(envir = .GlobalEnv)),
     envir = .GlobalEnv)
  suppressWarnings(source(configFile, local = FALSE))
  if (!exists("whichCompareVec") || is.null(whichCompareVec)) {
    cat(sprintf("[skip] %s: no whichCompareVec\n", sid))
    next
  }
  if (!exists("chanIntVec") || is.null(chanIntVec)) {
    cat(sprintf("[skip] %s: no chanIntVec\n", sid))
    next
  }

  dataFile <- file.path(dataDir, paste0(sid, "_PairedPulseData_avg_5.csv"))
  if (!file.exists(dataFile)) {
    cat(sprintf("[skip] %s: CSV missing\n", sid))
    next
  }
  dataPP <- read.table(dataFile, header = TRUE, sep = ",",
                       stringsAsFactors = FALSE,
                       colClasses = c("stimLevelVec" = "numeric",
                                      "sidVec" = "character"))
  dataPP$PPvec <- dataPP$PPvec * 1e6  # V -> uV
  dataPP <- subset(dataPP, chanVec %in% chanIntVec &
                           PPvec > 10 & PPvec < 1000)

  stimLevels <- sort(unique(dataPP$stimLevelVec))
  dataPP$mapStimLevel <- match(dataPP$stimLevelVec, stimLevels)

  for (ch in chanIntVec) {
    chData <- dataPP[dataPP$chanVec == ch, ]
    if (nrow(chData) == 0) next

    for (pair in whichCompareVec) {
      if (!is.numeric(pair) || length(pair) < 2) next
      pre_block  <- pair[1]
      post_block <- pair[2]

      post_idx <- which(blockIntLM == post_block)
      if (length(post_idx) == 0) next
      cond_raw <- blockType[post_idx]
      if (cond_raw == "baseline") next  # skip baseline-only pairs

      pairData <- chData[chData$blockVec %in% c(pre_block, post_block), ]
      if (nrow(pairData) == 0) next

      per_sl <- data.frame()
      for (sl in sort(unique(pairData$mapStimLevel))) {
        slData <- pairData[pairData$mapStimLevel == sl, ]
        pre_vals  <- slData$PPvec[slData$blockVec == pre_block]
        post_vals <- slData$PPvec[slData$blockVec == post_block]
        if (length(pre_vals) < 2 || length(post_vals) < 2) next

        obs_diff <- median(post_vals) - median(pre_vals)

        block_labels <- slData$blockVec
        all_vals <- slData$PPvec
        perm_diffs <- numeric(nPerm)
        for (p in 1:nPerm) {
          shuf <- sample(block_labels)
          perm_diffs[p] <- median(all_vals[shuf == post_block]) -
                            median(all_vals[shuf == pre_block])
        }
        p_val <- mean(abs(perm_diffs) >= abs(obs_diff))

        PPmed_pre  <- median(pre_vals)
        PPmed_post <- median(post_vals)
        pct_diff <- 100 * (PPmed_post - PPmed_pre) / PPmed_pre

        per_sl <- rbind(per_sl, data.frame(
          sid = sid,
          subjectNum = if (exists("subjectNum")) subjectNum else NA_integer_,
          chanVec = ch,
          pre_block = pre_block,
          post_block = post_block,
          condition = cond_raw,
          mapStimLevel = sl,
          stim_uA = stimLevels[sl],
          n_pre = length(pre_vals),
          n_post = length(post_vals),
          PPmed_pre_uV = round(PPmed_pre, 2),
          PPmed_post_uV = round(PPmed_post, 2),
          obs_diff_uV = round(obs_diff, 2),
          pct_diff = round(pct_diff, 1),
          perm_p = p_val,
          stringsAsFactors = FALSE
        ))
      }

      if (nrow(per_sl) > 0) {
        # BH-FDR within (subject, channel, condition) family across stim levels
        per_sl$perm_p_BH <- p.adjust(per_sl$perm_p, method = "BH")
        all_results[[length(all_results) + 1]] <- per_sl
        cat(sprintf("  %s ch %d %s (%d->%d): %d stim levels\n",
                    sid, ch, cond_raw, pre_block, post_block, nrow(per_sl)))
      }
    }
  }
}

perm_all <- do.call(rbind, all_results)

perm_all_out <- perm_all
perm_all_out$perm_p    <- signif(perm_all_out$perm_p, 4)
perm_all_out$perm_p_BH <- signif(perm_all_out$perm_p_BH, 4)

write.csv(perm_all_out,
          file.path(outDir, "permutation_results_all_subjects.csv"),
          row.names = FALSE)

perm_stim4 <- perm_all_out[perm_all_out$mapStimLevel == 4, ]
write.csv(perm_stim4,
          file.path(outDir, "permutation_results_all_subjects_stim4.csv"),
          row.names = FALSE)

cat("\n=== Summary (stim level 4 only) ===\n")
print(perm_stim4[order(perm_stim4$sid, perm_stim4$chanVec, perm_stim4$condition),
                 c("sid","chanVec","condition","pre_block","post_block",
                   "PPmed_pre_uV","PPmed_post_uV","pct_diff","perm_p","perm_p_BH")],
      row.names = FALSE)

cat(sprintf("\nSaved %d rows to permutation_results_all_subjects.csv\n", nrow(perm_all)))
cat(sprintf("Saved %d rows to permutation_results_all_subjects_stim4.csv\n", nrow(perm_stim4)))
