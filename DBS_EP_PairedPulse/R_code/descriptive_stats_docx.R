# Compute descriptive statistics for cell-median differences from baseline
# in raw uV (not logged). Saves to .docx for manuscript use.
# All values extracted from data — nothing hardcoded.

suppressPackageStartupMessages({
  library(plyr); library(here); library(dplyr); library(officer); library(flextable)
})

sidVec <- c('46c2a','c963f','2e114','fe7df','e6f3c','9f852',
            '8e907','08b13','e9c9b','41a73','68574','01fee','a23ed')

outputDir <- here("DBS_EP_PairedPulse","R_output")

allData <- list()
idx <- 1

for (sid in sidVec) {
  source(here("DBS_EP_PairedPulse","R_config_files",paste0("subj_",sid,".R")))
  dataPP <- read.table(here("DBS_EP_PairedPulse","R_data",paste0(sid,"_PairedPulseData_avg_5.csv")),
    header=TRUE, sep=",", stringsAsFactors=FALSE,
    colClasses=c("stimLevelVec"="numeric","sidVec"="factor"))
  dataPP$PPvec <- dataPP$PPvec * 1e6
  dataPP$stimLevelVec <- dataPP$stimLevelVec / 1e3
  dataPP <- subset(dataPP, chanVec %in% chanIntVec & PPvec > 10 & PPvec < 1000)
  dataPP <- na.exclude(dataPP)
  dataPP$blockVec <- as.factor(dataPP$blockVec)

  for (chanInt in chanIntVec) {
    dataInt <- subset(dataPP, chanVec == chanInt & blockVec %in% blockIntPlot)
    if (nrow(dataInt) == 0) next

    uniqueStimLevel <- sort(unique(dataInt$stimLevelVec))
    dataInt$mapStimLevel <- as.ordered(match(dataInt$stimLevelVec, uniqueStimLevel))
    dataInt$blockType <- mapvalues(dataInt$blockVec, from=blockIntPlot, to=blockType)

    for (comparison in whichCompareVec) {
      dc <- subset(dataInt, blockVec %in% comparison)
      dc <- dc %>%
        group_by(mapStimLevel, blockVec) %>% mutate(enough = n() >= 3) %>%
        group_by(mapStimLevel) %>% filter(all(enough) & length(unique(blockVec)) > 1) %>%
        ungroup()
      if (nrow(dc) == 0) next

      name <- as.character(unique(dc$blockType[dc$blockType != "baseline"]))
      if (length(name) == 0 || !name %in% c("A/A 200","A/A 25","A/B 200","A/B 25")) next

      # Cell medians in raw uV
      cellMed <- dc %>%
        group_by(mapStimLevel, blockVec, blockType) %>%
        summarise(PPmed = median(PPvec), .groups = "drop")

      bl <- cellMed %>% filter(blockType == "baseline") %>%
        dplyr::select(mapStimLevel, PPmed_bl = PPmed)
      post <- cellMed %>% filter(blockType != "baseline") %>%
        left_join(bl, by = "mapStimLevel") %>%
        mutate(
          diff_uV = PPmed - PPmed_bl,
          pct_diff = 100 * (PPmed - PPmed_bl) / PPmed_bl,
          condition = name, sid = sid, chan = chanInt)

      allData[[idx]] <- post
      idx <- idx + 1
    }
  }
}

df <- do.call(rbind, allData)

# --- Build summary tables ---
conditions <- c("A/A 200","A/A 25","A/B 200","A/B 25")

# Helper for summary stats
summarize_var <- function(x) {
  data.frame(
    n = length(x),
    Min = round(min(x), 1),
    Q25 = round(quantile(x, 0.25), 1),
    Median = round(median(x), 1),
    Q75 = round(quantile(x, 0.75), 1),
    Max = round(max(x), 1))
}

# Table 1: Absolute difference (uV) per cell
abs_rows <- lapply(conditions, function(cond) {
  d <- df %>% filter(condition == cond)
  cbind(Condition = cond, summarize_var(d$diff_uV))
})
abs_tbl <- do.call(rbind, abs_rows)

# Table 2: Percent difference per cell
pct_rows <- lapply(conditions, function(cond) {
  d <- df %>% filter(condition == cond)
  cbind(Condition = cond, summarize_var(d$pct_diff))
})
pct_tbl <- do.call(rbind, pct_rows)

# Table 3: Absolute difference, median across blocks per (subject, channel, stim level)
abs_block_rows <- lapply(conditions, function(cond) {
  d <- df %>% filter(condition == cond) %>%
    group_by(sid, chan, mapStimLevel) %>%
    summarise(val = median(diff_uV), .groups = "drop")
  cbind(Condition = cond, summarize_var(d$val))
})
abs_block_tbl <- do.call(rbind, abs_block_rows)

# Table 4: Percent difference, median across blocks
pct_block_rows <- lapply(conditions, function(cond) {
  d <- df %>% filter(condition == cond) %>%
    group_by(sid, chan, mapStimLevel) %>%
    summarise(val = median(pct_diff), .groups = "drop")
  cbind(Condition = cond, summarize_var(d$val))
})
pct_block_tbl <- do.call(rbind, pct_block_rows)

# --- Write .docx ---
doc <- read_docx()

doc <- body_add_par(doc, "Table: Absolute Difference from Baseline (uV, per cell median)", style = "heading 2")
ft <- flextable(abs_tbl) |> autofit() |>
  set_caption("Cell-median EP amplitude minus matched baseline cell median. Raw uV (not log-transformed).")
doc <- body_add_flextable(doc, ft)
doc <- body_add_par(doc, "")

doc <- body_add_par(doc, "Table: Percent Difference from Baseline (per cell median)", style = "heading 2")
ft <- flextable(pct_tbl) |> autofit() |>
  set_caption("100 * (post - baseline) / baseline, computed on cell-median raw uV values.")
doc <- body_add_flextable(doc, ft)
doc <- body_add_par(doc, "")

doc <- body_add_par(doc, "Table: Absolute Difference, Median Across Blocks (uV)", style = "heading 2")
ft <- flextable(abs_block_tbl) |> autofit() |>
  set_caption("Per (subject, channel, stim level), median of per-block absolute differences. Collapses repeated blocks within same condition.")
doc <- body_add_flextable(doc, ft)
doc <- body_add_par(doc, "")

doc <- body_add_par(doc, "Table: Percent Difference, Median Across Blocks", style = "heading 2")
ft <- flextable(pct_block_tbl) |> autofit() |>
  set_caption("Per (subject, channel, stim level), median of per-block percent differences.")
doc <- body_add_flextable(doc, ft)

outpath <- file.path(outputDir, "descriptive_stats_baseline_diff.docx")
print(doc, target = outpath)
cat("Saved to:", outpath, "\n")
