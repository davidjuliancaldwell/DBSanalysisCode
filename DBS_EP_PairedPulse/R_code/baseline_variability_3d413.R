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

min_stim_level = 2
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
    #dataPP <- subset(dataPP, PPvec<1000)
    #dataPP <- subset(dataPP, PPvec>30)
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
  labs(x = expression(paste("Stimulation Current (mA)")),y=expression(paste("log Peak to Peak Magnitude (",mu,"V)")),
       color="Experimental Condition",title = paste0("Anesthesia Effect on EP Magnitude"))
print(p3)

if (savePlot && !avgMeas) {
  ggsave(paste0("subj_", subjectNum, "_ID_", sid,"_combined.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps,fallback_resolution=600)
  ggsave(paste0("subj_", subjectNum, "_ID_", sid,"_combined.png"), units="in", width=figWidth, height=figHeight, dpi=600,device="png")
  
} else if (savePlot && avgMeas){
  ggsave(paste0("subj_", subjectNum, "_ID_", sid,"_combined_avg.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
  ggsave(paste0("subj_", subjectNum, "_ID_", sid,"_combined_avg.png"), units="in", width=figWidth, height=figHeight, dpi=600,device="png")
  
}

# fit models on trial-level data with block RE to handle pseudoreplication
# (1|blockVec) absorbs within-block correlation (10 blocks, sufficient for RE)
# chanVec as fixed effect (only 2 channels, too few for RE)

fit.lmm1 = lmerTest::lmer(PPvec ~ mapStimLevel + blockType + chanVec + (1|blockVec), data=dataList)
fit.lmm2 = lmerTest::lmer(PPvec ~ mapStimLevel + blockType * mapStimLevel + chanVec + (1|blockVec), data=dataList)

summary(fit.lmm1)
summary(glht(fit.lmm1,linfct=mcp(blockType="Tukey")))
summary(glht(fit.lmm1,linfct=mcp(mapStimLevel="Tukey")))

emmeans(fit.lmm1, list(pairwise ~ blockType), adjust = "tukey")
emmeans(fit.lmm1, list(pairwise ~ mapStimLevel), adjust = "tukey")

anova(fit.lmm1)
tab_model(fit.lmm1)

summary(fit.lmm2)
summary(glht(fit.lmm2,linfct=mcp(blockType="Tukey")))
summary(glht(fit.lmm2,linfct=mcp(mapStimLevel="Tukey")))

emmeans(fit.lmm2, list(pairwise ~ blockType), adjust = "tukey")
emmeans(fit.lmm2, list(pairwise ~ mapStimLevel), adjust = "tukey")

emm_s.t <- emmeans(fit.lmm2, pairwise ~ blockType | mapStimLevel)
emm_s.t <- emmeans(fit.lmm2, pairwise ~ mapStimLevel | blockType)

anova(fit.lmm2)
tab_model(fit.lmm2)

#########
# model comparison

anova(fit.lmm1, fit.lmm2)
lrtest(fit.lmm1, fit.lmm2)
AIC(fit.lmm1, fit.lmm2)
BIC(fit.lmm1, fit.lmm2)

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

setwd(oldwd)