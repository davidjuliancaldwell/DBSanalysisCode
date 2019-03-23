setwd('C:/Users/david/SharedCode/DBSanalysisCode')

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

rootDir = here()
dataDir = here("DBS_EP_PairedPulse","R_data")
codeDir = here("DBS_EP_PairedPulse","R_code")

sidVec = c("3d413")


savePlot = 0
avgMeasVec = c(0)
figWidth = 8 
figHeight = 4

dataList = list()
blockList = list()
conditionList = list()
index = 1

for (avgMeas in avgMeasVec) {
  for (sid in sidVec){
    source(here("DBS_EP_PairedPulse","R_config_files",paste0("subj_",sid,".R")))
    
    if (avgMeas) {
      dataPP <- read.table(here("DBS_EP_PairedPulse","R_data",paste0(sid,'_PairedPulseData_avg.csv')),header=TRUE,sep = ",",stringsAsFactors=F, colClasses=c("stimLevelVec"="numeric","sidVec"="character"))
    } else{
      dataPP <- read.table(here("DBS_EP_PairedPulse","R_data",paste0(sid,'_PairedPulseData.csv')),header=TRUE,sep = ",",stringsAsFactors=F, colClasses=c("stimLevelVec"="numeric","sidVec"="character"))
    }
    
    # multiply by 1e6
    dataPP$PPvec = dataPP$PPvec*1e6
    dataPP$stimLevelVec = dataPP$stimLevelVec/1e3
    dataPP <- subset(dataPP, PPvec<1500)
    dataPP <- subset(dataPP, PPvec>25)
    # change to factor 
    print(sid)
    
    for (chanInt in chanIntVec){
      
      blockType = blockTypeVec[[index]]
      blockIntPlot = blockIntPlotVec[[index]]
      # select data of interest 
      dataInt <- na.exclude(subset(dataPP,(chanVec == chanInt) & (blockVec %in% blockIntPlot)))
      dataInt$blockVec = as.factor(dataInt$blockVec)
      dataInt$chanVec = as.factor(dataInt$chanVec)
      
      # map stimulation levels to consistent ordering for between subject comparison
      uniqueStimLevel = as.double(unique(dataInt$stimLevelVec))
      
      
      mappingStimLevel =c(1:length(uniqueStimLevel))
      
      dataInt$mapStimLevel <- mapvalues(dataInt$stimLevelVec,
                                        from=uniqueStimLevel,
                                        to=mappingStimLevel)
      uniqueBlockLevel = unique(dataInt$blockVec)
      
      dataInt$blockType <- mapvalues(dataInt$blockVec,
                                     from=uniqueBlockLevel,
                                     to=blockType)
      
      dataInt$mapStimLevel = as.ordered(dataInt$mapStimLevel)
      dataInt$blockType = as.factor(dataInt$blockType)
      dataIntCompare <- dataInt
      
      dataIntCompare = as_data_frame(dataIntCompare)
      dataIntCompare$index = index
      dataList[[index]] = dataIntCompare

      fit.lm    = lm(PPvec ~ mapStimLevel + blockVec + mapStimLevel*blockVec,data=dataIntCompare)
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
  labs(x = expression(paste("Stimulation Current (mA)")),y=expression(paste("Peak to peak magnitude (",mu,"V)")),
       color="Experimental Condition",title = paste0("Anesthesia Effect on EP Magnitude"))
print(p3)

if (savePlot && !avgMeas) {
  ggsave(paste0("subj_", subjectNum, "_ID_", sid,"_combined.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps,fallback_resolution=600)
  ggsave(paste0("subj_", subjectNum, "_ID_", sid,"_combined.png"), units="in", width=figWidth, height=figHeight, dpi=600,device="png")
  
} else if (savePlot && avgMeas){
  ggsave(paste0("subj_", subjectNum, "_ID_", sid,"_combined_avg.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
  ggsave(paste0("subj_", subjectNum, "_ID_", sid,"_combined_avg.png"), units="in", width=figWidth, height=figHeight, dpi=600,device="png")
  
}

# fit models

fit.lmm = lmerTest::lmer(PPvec ~ mapStimLevel + blockType + (1|chanVec),data=dataList)

fit.lmm2 = lmerTest::lmer(PPvec ~ mapStimLevel + blockType + blockType*mapStimLevel + (1|chanVec),data=dataList)

summary(fit.lmm)
# plot(fit.lm)
summary(glht(fit.lmm,linfct=mcp(blockType="Tukey")))
summary(glht(fit.lmm,linfct=mcp(mapStimLevel="Tukey")))

emmeans(fit.lmm, list(pairwise ~ blockType), adjust = "tukey")
emmeans(fit.lmm, list(pairwise ~ mapStimLevel), adjust = "tukey")

#emm_s.t <- emmeans(fit.lmm, pairwise ~ blockType | mapStimLevel)
#emm_s.t <- emmeans(fit.lmm, pairwise ~ mapStimLevel | blocktype)

anova(fit.lmm)
tab_model(fit.lmm)

summary(fit.lmm2)
# plot(fit.lm)
summary(glht(fit.lmm2,linfct=mcp(blockType="Tukey")))
summary(glht(fit.lmm2,linfct=mcp(mapStimLevel="Tukey")))

emmeans(fit.lmm2, list(pairwise ~ blockType), adjust = "tukey")
emmeans(fit.lmm2, list(pairwise ~ mapStimLevel), adjust = "tukey")

emm_s.t <- emmeans(fit.lmm2, pairwise ~ blockType | mapStimLevel)
emm_s.t <- emmeans(fit.lmm2, pairwise ~ mapStimLevel | blocktype)

anova(fit.lmm2)
tab_model(fit.lmm2)



#########
# make some plots!

anova(fit.lmm, fit.lmm2)
# Likelihood ratio test

lrtest(fit.lmm, fit.lmm2)

# aic
AIC(fit.lmm, fit.lmm2)

BIC(fit.lmm, fit.lmm2)