library('plyr')
library('here')
library('nlme')
library('ggplot2')
library('drc')
library('minpack.lm')
library('lmtest')
library('glmm')
library("lme4")
library('multcomp')
library('lmerTest')
library('sjPlot')
library('emmeans')
library("dplyr")
library("DescTools")
library('effectsize')

rootDir = here()
dataDir = here("DBS_EP_PairedPulse","R_data")
codeDir = here("DBS_EP_PairedPulse","R_code")
outputDir = here("DBS_EP_PairedPulse","R_output")
dir.create(outputDir, showWarnings = FALSE)

sidVec = c("3d413")

min_stim_level = 1
log_data = FALSE
box_data = FALSE
trim_data = FALSE
savePlot = 1
avgMeasVec = c(1)
figWidth = 8
figHeight = 4

dataList = list()
blockList = list()
conditionList = list()
index = 1

oldwd <- getwd()
setwd(outputDir)

for (avgMeas in avgMeasVec) {
  for (sid in sidVec){
    source(here("DBS_EP_PairedPulse","R_config_files",paste0("subj_",sid,".R")))
    brodmann_areas <- read.csv(here("DBS_EP_PairedPulse","R_config_files",paste0(sid,'_MNIcoords_labelled.csv')),header=TRUE,sep = ",",stringsAsFactors=F)
    
    dataPP <- read.table(here("DBS_EP_PairedPulse","R_data",paste0(sid,'_PairedPulseData_avg_5.csv')),header=TRUE,sep = ",",stringsAsFactors=F, colClasses=c("stimLevelVec"="numeric","sidVec"="character"))

    # PPvec is in volts, convert to microvolts
    dataPP$PPvec = dataPP$PPvec*1e6
    dataPP$stimLevelVec = dataPP$stimLevelVec/1e3
    dataPP <- subset(dataPP, PPvec<1000)
    dataPP <- subset(dataPP, PPvec>30)
    # change to factor 
    print(sid)
    
    for (chanInt in chanIntVec){
      
      blockType = blockTypeVec[[index]]
      blockIntPlot = blockIntPlotVec[[index]]
      # select data of interest 
      dataInt <- subset(dataPP,(chanVec == chanInt) & (blockVec %in% blockIntPlot))
      dataInt$blockVec = as.factor(dataInt$blockVec)
      dataInt$chanVec = as.factor(dataInt$chanVec)
      
      # map stimulation levels to consistent ordering for between subject comparison
      uniqueStimLevel = as.double(unique(dataInt$stimLevelVec))
      
      if (sid == "2e114"){
        mappingStimLevel =c(1,3,4)
        
      } else {
        mappingStimLevel =c(1:length(uniqueStimLevel))
      }
      dataInt$mapStimLevel<- mapvalues(dataInt$stimLevelVec,
                                       from=sort(uniqueStimLevel),
                                       to=mappingStimLevel)
      
     # uniqueBlockLevel = unique(dataInt$blockVec)
     # blockTypeTrim = blockType[uniqueBlockLevel]
      
     # dataInt$blockType <- mapvalues(dataInt$blockVec,
     #                                from=uniqueBlockLevel,
     #                                to=blockTypeTrim)
      
      
      dataInt$blockType <- mapvalues(dataInt$blockVec,
                                     from=blockIntPlot,
                                     to=blockType)
      
      dataInt$mapStimLevel = as.ordered(dataInt$mapStimLevel)
      dataInt$blockType = as.factor(dataInt$blockType)
      
      
      dataInt <- na.exclude(dataInt)
      
      if (trim_data){
        dataInt <- dataInt %>% group_by(blockVec,blockType,mapStimLevel) %>% mutate(PPvecLabel = !is.element(seq_len(length(PPvec)),attr(Trim(PPvec,0.025,na.rm=TRUE),'trim')))
        dataInt <- subset(dataInt,PPvecLabel==TRUE)
        
      }
      
      dataInt <- subset(dataInt,mapStimLevel>=min_stim_level)
      
      if (log_data){
        dataInt$PPvec <- log(dataInt$PPvec)
      }
      
      if (box_data){
        box_trans <- caret::BoxCoxTrans(dataInt$PPvec)
        dataInt$PPvec <- predict(box_trans,dataInt$PPvec)
      }
      
      # get which brodmann area label is for this electrode
      dataInt$baLabel = as.character(brodmann_areas$ba.label[chanInt])
      dataInt$aalLabel = as.character(brodmann_areas$aal.label[chanInt])
      
      # get all to same side
      dataInt$baLabel[dataInt$baLabel=='Right-PrimMotor (4)'] = 'Left-PrimMotor (4)'
      dataInt$aalLabel[dataInt$aalLabel=='Postcentral_R'] = 'Postcentral_L'
      
      dataInt$baLabel = as.factor(dataInt$baLabel)
      dataInt$aalLabel = as.factor(dataInt$aalLabel)
      
      dataIntCompare <- dataInt
      
      dataIntCompare = as_data_frame(dataIntCompare)
      dataIntCompare$index = index
      
      dataList[[index]] = dataIntCompare

      #fit.lm    = lm(PPvec ~ mapStimLevel + blockVec + mapStimLevel*blockVec,data=dataIntCompare)
      fit.lm    = lm(PPvec ~ mapStimLevel + blockVec,data=dataIntCompare)
      
      summary(fit.lm)
      # plot(fit.lm)
      summary(glht(fit.lm,linfct=mcp(blockVec="Tukey")))
      summary(glht(fit.lm,linfct=mcp(mapStimLevel="Tukey")))
      
      emmeans(fit.lm, list(pairwise ~ blockVec), adjust = "tukey")
      emmeans(fit.lm, list(pairwise ~ mapStimLevel), adjust = "tukey")
      
      emm_s.t <- emmeans(fit.lm, pairwise ~ blockVec | mapStimLevel)
      emm_s.t <- emmeans(fit.lm, pairwise ~ mapStimLevel | blockVec)
      
      anova(fit.lm)
      tab_model(fit.lm)
      
      p <- ggplot(dataInt, aes(x=stimLevelVec, y=PPvec,color=stimLevelVec)) +
        geom_point(position=position_jitterdodge(dodge.width=0.250)) +  geom_smooth(method=lm) + facet_wrap(~blockVec, labeller = as_labeller(blockNames))+
        labs(x = expression(paste("Stimulation Current (mA)")),y=expression(paste("Peak to Peak magnitude (",mu,"V)")),title = paste0("Subject ", subjectNum, " ID ", sid," DBS Paired Pulse EP Mseasurements")) +
        guides(colour=guide_colorbar("Stimulation Level"))
      print(p)
      
      if (savePlot && !avgMeas) {
        ggsave(paste0("subj_", subjectNum, "_ID_", sid,"_scatter_lm.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
         ggsave(paste0("subj_", subjectNum, "_ID_", sid,"_scatter_lm.png"), units="in", width=figWidth, height=figHeight, dpi=600)
        
      }
      else if (savePlot && avgMeas){
        ggsave(paste0("subj_", subjectNum, "_ID_", sid,"_scatter_lm_avg.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
        ggsave(paste0("subj_", subjectNum, "_ID_", sid,"_scatter_lm_avg.png"), units="in", width=figWidth, height=figHeight, dpi=600)
        
      }
      index = index + 1 
    }
  }
}
############# 
# now do comparison at highest stim level relative to baseline
dataList = do.call(rbind,dataList)

#plot

labelVec <- c("4" = "Channel 4","6" = "Channel 6")

#dataListFactorized$stimLevelVec = as.factor(dataListFactorized$stimLevelVec)
p3 <- ggplot(dataList, aes(x=stimLevelVec, y=PPvec,color=blockType)) +
  geom_point(position=position_jitterdodge(dodge.width=0.250)) +geom_smooth(method=lm) + facet_grid(. ~chanVec, labeller = labeller(chanVec = labelVec))+
  labs(x = expression(paste("Stimulation Current (mA)")),y=expression(paste("Peak to Peak Magnitude (",mu,"V)")),
       color="Experimental Condition",title = paste0("Anesthesia Effect on EP Magnitude"))
print(p3)

if (savePlot && !avgMeas) {
  ggsave(paste0("subj_", subjectNum, "_ID_", sid,"_combined.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps,fallback_resolution=600)
  ggsave(paste0("subj_", subjectNum, "_ID_", sid,"_combined.png"), units="in", width=figWidth, height=figHeight, dpi=600,device="png")
  
} else if (savePlot && avgMeas){
  ggsave(paste0("subj_", subjectNum, "_ID_", sid,"_combined_avg.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
  ggsave(paste0("subj_", subjectNum, "_ID_", sid,"_combined_avg.png"), units="in", width=figWidth, height=figHeight, dpi=600,device="png")
  
}

cat(sprintf("Trial-level observations: %d\n", nrow(dataList)))

# Trial-level lmer with uncorrelated random linear dose slope per block.
# The random slope allows each block to have its own dose-response
# steepness, fixing stim-level pseudoreplication (df drops from ~459
# to ~23) while preserving correct ~7 df for awake/asleep.
# stim_linpoly is the linear polynomial contrast from contr.poly().
n_stim_levels <- nlevels(dataList$mapStimLevel)
poly_lin <- contr.poly(n_stim_levels)[, 1]
dataList$stim_linpoly <- poly_lin[as.numeric(dataList$mapStimLevel)]

# Sum-to-zero coding for blockType (used in interaction model fit.lmm2).
# mapStimLevel is ordered (contr.poly, already centered).
contrasts(dataList$blockType) <- contr.sum

fit.lmm1 = lmerTest::lmer(PPvec ~ mapStimLevel + blockType + chanVec +
  (1|blockVec) + (0+stim_linpoly|blockVec), data=dataList)
fit.lmm2 = lmerTest::lmer(PPvec ~ mapStimLevel * blockType + chanVec +
  (1|blockVec) + (0+stim_linpoly|blockVec), data=dataList)

cat("\n=== Primary model: additive with random slope ===\n")
print(summary(fit.lmm1))
cat("\n=== ANOVA (additive) ===\n")
print(anova(fit.lmm1, type=3))

emmeans(fit.lmm1, list(pairwise ~ blockType), adjust = "tukey")
emmeans(fit.lmm1, list(pairwise ~ mapStimLevel), adjust = "tukey")

tab_model(fit.lmm1)

cat("\n=== Interaction model ===\n")
print(summary(fit.lmm2))
cat("\n=== ANOVA (interaction) ===\n")
print(anova(fit.lmm2, type=3))

tab_model(fit.lmm2)

cat("\n=== Model comparison: additive vs interaction ===\n")
print(anova(fit.lmm1, fit.lmm2, refit=FALSE))
print(AIC(fit.lmm1, fit.lmm2))

#########
# Effect size analysis (appropriate for single-subject design)
# Cohen's d per (channel, stim_level): awake vs asleep
# NOTE: PPvec is log-transformed (log_data=TRUE), so d is on log scale

es_list <- list()
es_idx <- 1
for (chan in unique(dataList$chanVec)) {
  chanData <- dataList %>% filter(chanVec == chan)
  for (sl in unique(chanData$mapStimLevel)) {
    slData <- chanData %>% filter(mapStimLevel == sl)
    awake <- slData$PPvec[slData$blockType == "awake"]
    asleep <- slData$PPvec[slData$blockType == "asleep"]
    if (length(awake) >= 5 & length(asleep) >= 5) {
      d <- cohens_d(awake, asleep)
      es_list[[es_idx]] <- data.frame(
        chanVec = chan,
        mapStimLevel = sl,
        effect.size = d$Cohens_d,
        CI_low = d$CI_low,
        CI_high = d$CI_high
      )
      es_idx <- es_idx + 1
    }
  }
}
es_df <- do.call(rbind, es_list)

# --- Effect size bar plot: awake vs asleep by stim level and channel ---
p_es <- ggplot(es_df, aes(x = mapStimLevel, y = effect.size, fill = chanVec)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high),
                position = position_dodge(width = 0.7), width = 0.2) +
  geom_hline(yintercept = c(-0.8, -0.5, -0.2, 0, 0.2, 0.5, 0.8),
             linetype = "dashed", color = "grey50", linewidth = 0.3) +
  labs(x = "Stimulation Level",
       y = "Cohen's d (95% CI)",
       fill = "Channel",
       title = "Effect Size: Awake vs Asleep by Stimulation Level",
       subtitle = "Computed on log-transformed EP magnitude") +
  theme_bw(base_size = 14)
print(p_es)

if (savePlot) {
  ggsave(paste0("subj_3d413_effect_size_awake_asleep.png"),
         plot = p_es, units = "in", width = figWidth, height = figHeight, dpi = 600)
  ggsave(paste0("subj_3d413_effect_size_awake_asleep.eps"),
         plot = p_es, units = "in", width = figWidth, height = figHeight,
         device = cairo_ps, fallback_resolution = 600)
}

#########
# Permutation tests (trial-level, per channel)
# Shuffles awake/asleep labels across trials within each stim level
# Stratified pooled test: independently permute within each stim level,
#   then average per-stim-level median differences
# Consistent with a23ed approach; note slight liberality due to ICC ~ 0.07
# Uses median as test statistic

set.seed(42)
nPerm <- 10000

# --- Per channel, per stim level ---
perm_results <- data.frame(channel = character(), mapStimLevel = character(),
                           obs_diff = numeric(), perm_p = numeric(),
                           stringsAsFactors = FALSE)

for (chan in unique(dataList$chanVec)) {
  chanData <- dataList %>% filter(chanVec == chan)
  for (sl in unique(chanData$mapStimLevel)) {
    slData <- chanData %>% filter(mapStimLevel == sl)
    obs_diff <- median(slData$PPvec[slData$blockType == "awake"]) -
                median(slData$PPvec[slData$blockType == "asleep"])

    labels <- slData$blockType
    vals <- slData$PPvec
    perm_diffs <- numeric(nPerm)
    for (p in 1:nPerm) {
      shuf <- sample(labels)
      perm_diffs[p] <- median(vals[shuf == "awake"]) - median(vals[shuf == "asleep"])
    }
    p_val <- mean(abs(perm_diffs) >= abs(obs_diff))
    perm_results <- rbind(perm_results,
      data.frame(channel = as.character(chan), mapStimLevel = as.character(sl),
                 obs_diff = obs_diff, perm_p = p_val))
  }
}

cat("\n=== Permutation tests (median): awake vs asleep (per channel, per stim level) ===\n")
print(perm_results)

# --- Per channel, stratified pooled across stim levels ---
# Independently permute within each stim level, then average the per-stim-level
# median differences. Preserves stim level balance in every permutation.
perm_stratified <- data.frame(channel = character(), obs_diff = numeric(),
                              perm_p = numeric(), stringsAsFactors = FALSE)
chan_perm_null <- list()

for (chan in unique(dataList$chanVec)) {
  chanData <- dataList %>% filter(chanVec == chan)
  stim_levels <- sort(unique(chanData$mapStimLevel))

  # Observed stratified statistic
  obs_sl_diffs <- numeric(length(stim_levels))
  for (j in seq_along(stim_levels)) {
    slData <- chanData %>% filter(mapStimLevel == stim_levels[j])
    obs_sl_diffs[j] <- median(slData$PPvec[slData$blockType == "awake"]) -
                       median(slData$PPvec[slData$blockType == "asleep"])
  }
  obs_pooled <- mean(obs_sl_diffs)

  # Permutation null: independently shuffle within each stim level
  perm_pooled <- numeric(nPerm)
  for (p in 1:nPerm) {
    perm_sl_diffs <- numeric(length(stim_levels))
    for (j in seq_along(stim_levels)) {
      slData <- chanData %>% filter(mapStimLevel == stim_levels[j])
      shuf <- sample(slData$blockType)
      perm_sl_diffs[j] <- median(slData$PPvec[shuf == "awake"]) -
                          median(slData$PPvec[shuf == "asleep"])
    }
    perm_pooled[p] <- mean(perm_sl_diffs)
  }
  pooled_p <- mean(abs(perm_pooled) >= abs(obs_pooled))

  perm_stratified <- rbind(perm_stratified,
    data.frame(channel = as.character(chan), obs_diff = obs_pooled, perm_p = pooled_p))
  chan_perm_null[[as.character(chan)]] <- perm_pooled

  cat(sprintf("\nChannel %s stratified pooled: obs = %.2f uV, p = %.4f\n",
              chan, obs_pooled, pooled_p))
}

cat("\n=== Permutation tests (median): awake vs asleep (stratified pooled) ===\n")
print(perm_stratified)

# --- Null distribution histograms (per channel, stratified pooled) ---
for (chan in unique(dataList$chanVec)) {
  obs_val <- perm_stratified$obs_diff[perm_stratified$channel == as.character(chan)]
  p_val <- perm_stratified$perm_p[perm_stratified$channel == as.character(chan)]
  null_df <- data.frame(value = chan_perm_null[[as.character(chan)]])

  p_perm <- ggplot(null_df, aes(x = value)) +
    geom_histogram(aes(fill = abs(value) >= abs(obs_val)),
                   bins = 60, color = "grey50", show.legend = FALSE) +
    scale_fill_manual(values = c("FALSE" = "grey70", "TRUE" = "#C93312")) +
    geom_vline(xintercept = obs_val, color = "#C93312", linewidth = 1.2) +
    geom_vline(xintercept = -obs_val, color = "#C93312",
               linewidth = 1.2, linetype = "dashed") +
    annotate("text", x = obs_val, y = Inf, vjust = 2, hjust = -0.1,
             label = sprintf("observed = %.1f uV\np = %.4f", obs_val, p_val),
             color = "#C93312", size = 4) +
    labs(x = expression(paste("Awake - Asleep  [", mu, "V]")),
         y = "Count",
         title = sprintf("Permutation Test: Awake vs Asleep (Channel %s, Stratified)", chan),
         subtitle = sprintf("%s permutations, stratified by stim level (two-sided)",
                            formatC(nPerm, format = "d", big.mark = ","))) +
    theme_bw(base_size = 14)
  print(p_perm)

  if (savePlot) {
    ggsave(sprintf("subj_3d413_perm_awake_asleep_chan%s.png", chan),
           plot = p_perm, units = "in", width = 7, height = 4, dpi = 600)
    ggsave(sprintf("subj_3d413_perm_awake_asleep_chan%s.eps", chan),
           plot = p_perm, units = "in", width = 7, height = 4,
           device = cairo_ps, fallback_resolution = 600)
  }
}

# --- Save manuscript-ready tables ---
# tab_model HTML (Word-compatible)
tab_3d413 <- tab_model(fit.lmm1, fit.lmm2, show.re.var=TRUE, show.icc=TRUE, show.obs=TRUE,
  dv.labels=c("Additive (primary)", "Interaction"),
  title="3d413: Awake vs Asleep (trial-level, random dose slope)")
writeLines(as.character(tab_3d413$page.complete), paste0(outputDir, "/tab_model_3d413.html"))
cat("Saved tab_model HTML to:", paste0(outputDir, "/tab_model_3d413.html"), "\n")

# .docx with ANOVA, fixed effects, awake/asleep contrast, and permutation results
if (requireNamespace("officer", quietly=TRUE) && requireNamespace("flextable", quietly=TRUE)) {
  library(officer)
  library(flextable)

  fmt_p <- function(p) ifelse(p < 0.001, "< 0.001", sprintf("%.3f", p))
  doc <- read_docx()

  # ANOVA (lmerTest Type III with Satterthwaite df)
  aov_tbl <- as.data.frame(anova(fit.lmm1, type=3))
  aov_tbl$Effect <- rownames(aov_tbl)
  aov_tbl <- aov_tbl[, c("Effect","Sum Sq","Mean Sq","NumDF","DenDF","F value","Pr(>F)")]
  aov_tbl$`Sum Sq` <- round(aov_tbl$`Sum Sq`, 1)
  aov_tbl$`Mean Sq` <- round(aov_tbl$`Mean Sq`, 1)
  aov_tbl$DenDF <- round(aov_tbl$DenDF, 1)
  aov_tbl$`F value` <- round(aov_tbl$`F value`, 2)
  aov_tbl$p <- sapply(aov_tbl$`Pr(>F)`, fmt_p)
  aov_tbl$`Pr(>F)` <- NULL
  doc <- body_add_par(doc, "Table: ANOVA (3d413, additive model, random slope)", style="heading 2")
  ft <- flextable(aov_tbl) |> autofit() |>
    set_caption("Trial-level lmer with uncorrelated random linear dose slope per block")
  doc <- body_add_flextable(doc, ft)
  doc <- body_add_par(doc, "")

  # Awake vs asleep emmeans
  doc <- body_add_par(doc, "Table: Awake vs Asleep (3d413)", style="heading 2")
  emm_bt <- emmeans(fit.lmm1, pairwise ~ blockType)
  bt_tbl <- as.data.frame(emm_bt$contrasts)
  bt_tbl$estimate <- round(bt_tbl$estimate, 2)
  bt_tbl$SE <- round(bt_tbl$SE, 2)
  bt_tbl$df <- round(bt_tbl$df, 1)
  bt_tbl$t.ratio <- round(bt_tbl$t.ratio, 3)
  bt_tbl$p.value <- sapply(bt_tbl$p.value, fmt_p)
  ft <- flextable(bt_tbl) |> autofit() |>
    set_caption("Awake vs asleep emmeans contrast from additive model")
  doc <- body_add_flextable(doc, ft)
  doc <- body_add_par(doc, "")
  # Random effects
  doc <- body_add_par(doc, "Table: Random Effects (3d413)", style="heading 2")
  vc <- VarCorr(fit.lmm1)
  ngrps <- summary(fit.lmm1)$ngrps
  re_rows <- list()
  for (nm in names(vc)) {
    v <- vc[[nm]]
    ng <- as.character(ngrps[nm])
    re_rows[[length(re_rows)+1]] <- data.frame(
      Component=nm, Term=rownames(as.matrix(v)),
      Variance=round(diag(as.matrix(v)),2), SD=round(sqrt(diag(as.matrix(v))),2),
      Groups=c(ng, rep("", ncol(as.matrix(v))-1)), stringsAsFactors=FALSE)
  }
  re_df <- do.call(rbind, re_rows)
  re_df <- rbind(re_df, data.frame(Component="Residual", Term="",
    Variance=round(sigma(fit.lmm1)^2,2), SD=round(sigma(fit.lmm1),2), Groups=""))
  ft <- flextable(re_df) |> autofit() |>
    set_caption("Random effects: block intercept + uncorrelated linear dose slope per block")
  doc <- body_add_flextable(doc, ft)
  doc <- body_add_par(doc, "")

  # Fixed effects
  fe_tbl <- as.data.frame(summary(fit.lmm1)$coefficients)
  fe_tbl$Predictor <- rownames(fe_tbl)
  fe_tbl <- fe_tbl[, c("Predictor","Estimate","Std. Error","df","t value","Pr(>|t|)")]
  fe_tbl$Estimate <- round(fe_tbl$Estimate, 2)
  fe_tbl$`Std. Error` <- round(fe_tbl$`Std. Error`, 2)
  fe_tbl$df <- round(fe_tbl$df, 1)
  fe_tbl$`t value` <- round(fe_tbl$`t value`, 2)
  fe_tbl$p <- sapply(fe_tbl$`Pr(>|t|)`, fmt_p)
  fe_tbl$`Pr(>|t|)` <- NULL
  doc <- body_add_par(doc, "Table: Fixed Effects (3d413, additive model)", style="heading 2")
  ft <- flextable(fe_tbl) |> autofit()
  doc <- body_add_flextable(doc, ft)
  doc <- body_add_par(doc, "Raw uV scale. Trial-level lmer with random dose slope per block.")
  doc <- body_add_par(doc, "")

  # Permutation results
  doc <- body_add_par(doc, "Table: Permutation Tests (3d413)", style="heading 2")
  perm_results$obs_diff <- round(perm_results$obs_diff, 2)
  perm_results$perm_p <- sapply(perm_results$perm_p, fmt_p)
  ft <- flextable(perm_results) |> autofit() |>
    set_caption("Per channel, per stim level (10,000 permutations, median statistic)")
  doc <- body_add_flextable(doc, ft)
  doc <- body_add_par(doc, "")

  doc <- body_add_par(doc, "Table: Stratified Pooled Permutation Tests (3d413)", style="heading 2")
  perm_stratified$obs_diff <- round(perm_stratified$obs_diff, 2)
  perm_stratified$perm_p <- sapply(perm_stratified$perm_p, fmt_p)
  ft <- flextable(perm_stratified) |> autofit() |>
    set_caption("Stratified pooled per channel (10,000 permutations)")
  doc <- body_add_flextable(doc, ft)

  docx_path <- paste0(outputDir, "/statistical_tables_3d413.docx")
  print(doc, target=docx_path)
  cat("Saved 3d413 manuscript tables to:", docx_path, "\n")
}

setwd(oldwd)