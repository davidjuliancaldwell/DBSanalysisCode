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

sidVec <- c('46c2a','c963f','2e114','9f852',
            '08b13','8e907','e9c9b','41a73','68574',
            '01fee','a23ed')

# sidVec <- c('46c2a','c963f','2e114',
#             '08b13','8e907','e9c9b','41a73','68574',
#             '01fee')

#sidVec <- c('9f852')
#sidVec <- c('46c2a')

savePlot = 0
avgMeasVec = c(0)
figWidth = 8 
figHeight = 6 


for (avgMeas in avgMeasVec) {
  
  dataList = list()
  blockList = list()
  conditionList = list()
  index = 1
  
  for (sid in sidVec){
    source(here("DBS_EP_PairedPulse","R_config_files",paste0("subj_",sid,".R")))
    
    if (avgMeas) {
      dataPP <- read.table(here("DBS_EP_PairedPulse","R_data",paste0(sid,'_PairedPulseData_avg.csv')),header=TRUE,sep = ",",stringsAsFactors=F, colClasses=c("stimLevelVec"="numeric","sidVec"="character"))
    } else{
      dataPP <- read.table(here("DBS_EP_PairedPulse","R_data",paste0(sid,'_PairedPulseData.csv')),header=TRUE,sep = ",",stringsAsFactors=F, colClasses=c("stimLevelVec"="numeric","sidVec"="character"))
    }
    
    if (sid == "2e114"){
      goodVec = c(1000,2500,3500)
      dataPP <- dataPP %>% filter(stimLevelVec %in% goodVec)
    }
    
    # multiply by 1e6
    dataPP$PPvec = dataPP$PPvec*1e6
    dataPP$stimLevelVec = dataPP$stimLevelVec/1e3
    dataPP <- subset(dataPP, PPvec<1500)
    dataPP <- subset(dataPP, PPvec>25)
    
    # denote which channel was in the conditioning pair 
    
    dataPP <- subset(dataPP,(chanVec %in% chanIntVec))
    
    dataPP$chanInCond <- mapvalues(dataPP$chanVec,
                                   from=chanIntVec,
                                   to=chanIntConditioningPairVec)
    # change to factor 
    dataPP$blockVec = as.factor(dataPP$blockVec)
    print(sid)
    
    for (chanInt in chanIntVec){
      
      # select data of interest 
      dataInt <- na.exclude(subset(dataPP,(chanVec == chanInt) & (blockVec %in% blockIntPlot)))
      
      # map stimulation levels to consistent ordering for between subject comparison
      uniqueStimLevel = as.double(unique(dataInt$stimLevelVec))
      
      if (sid == "2e114"){
        mappingStimLevel =ordered(c(1,3,4))
        
      } else {
        mappingStimLevel =ordered(c(1:length(uniqueStimLevel)))
      }
      dataInt$mapStimLevel <- mapvalues(dataInt$stimLevelVec,
                                        from=uniqueStimLevel,
                                        to=mappingStimLevel)
      
      uniqueBlockLevel = unique(dataInt$blockVec)
      blockTypeTrim = blockType[uniqueBlockLevel]
      
      dataInt$blockType <- mapvalues(dataInt$blockVec,
                                     from=uniqueBlockLevel,
                                     to=blockTypeTrim)
      
      dataInt$mapStimLevel = as.factor(dataInt$mapStimLevel)
      dataInt$blockType = as.factor(dataInt$blockType)
      for (comparison in whichCompareVec){
        dataIntCompare <- subset(dataInt,(blockVec %in% comparison))
        
        dataIntCompare <- dataIntCompare %>% 
          group_by(mapStimLevel) %>% 
          mutate(absDiff = PPvec - mean(PPvec[blockVec==comparison[1]]),
                 percentDiff = 100*(PPvec - mean(PPvec[blockVec==comparison[1]]))/mean(PPvec[blockVec==comparison[1]]) ) 
        
        dataIntCompare = as_data_frame(dataIntCompare)
        dataIntCompare$index = index
        dataList[[index]] = dataIntCompare
        blockList[[index]] = comparison
        
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
        
        index = index + 1
        
      }
      
      p <- ggplot(dataInt, aes(x=stimLevelVec, y=PPvec,colour=stimLevelVec)) +
        geom_point(position=position_jitterdodge(dodge.width=0.250)) +  geom_smooth(method=lm) + facet_wrap(~blockVec, labeller = as_labeller(blockNames))+
        labs(x = expression(paste("Stimulation current (mA)")),y=expression(paste("Peak to Peak magnitude (",mu,"V)"), fill="stimulus level"),title = paste0("Subject ", subjectNum, " ID ", sid," DBS Paired Pulse EP Measurements")) +
        guides(colour=guide_colorbar("stimulation level"))
      print(p)
      
      if (savePlot && !avgMeas) {
        ggsave(paste0("subj_", subjectNum, "_ID_", sid,"_chan_",chanInt,"_scatter_lm.png"), units="in", width=figWidth, height=figHeight, dpi=600)
        ggsave(paste0("subj_", subjectNum, "_ID_", sid,"_chan_",chanInt,"_scatter_lm.eps"), units="in", width=figWidth, height=figHeight, dpi=600,device="eps")
        
      }
      else if (savePlot && avgMeas){
        ggsave(paste0("subj_", subjectNum, "_ID_", sid,"_chan_",chanInt,"_scatter_lm_avg.png"), units="in", width=figWidth, height=figHeight, dpi=600)
        ggsave(paste0("subj_", subjectNum, "_ID_", sid,"_chan_",chanInt,"_scatter_lm_avg.eps"), units="in", width=figWidth, height=figHeight, dpi=600,device="eps")
        
      }
      
      dataToFit <- na.exclude(dataPP[dataPP$blockVec %in% blockIntLM & dataPP$chanVec==chanInt,])
      fit.glm    = glm(PPvec ~ blockVec + stimLevelVec,data=dataToFit)
      summary(fit.glm)
      summary(glht(fit.glm,linfct=mcp(blockVec="Tukey")))
      # plot(fit.glm)
      # 
      # # fit.lmm = lmer(PPvec ~ stimLevelVec + blockVec + (1|stimLevelVec) + (1|blockVec), data=dataToFit)
      # # summary(fit.lmm)
      # #confint(fit.lmm,method="boot")
      # 
      # # get starting values
      # if (sid != '8e907') {
      #   fit.startnlme0 <- nlsLM(PPvec ~ SSlogis(stimLevelVec, Asym, xmid, scal),data=dataToFit)
      # } else {
      #   fit.startnlme0 <- nlsLM(PPvec ~ SSlogis(stimLevelVec, Asym, xmid, scal),data=dataToFit,start = c(Asym=400,xmid = 500,scal = 500),control = list(maxiter = 500))
      # }
      # 
      # startVals = summary(fit.startnlme0)$coefficients
      # # 
      # #  if (sid == '08b13') {
      # #    fit.nlme0 <- gnls(PPvec ~ SSlogis(stimLevelVec, Asym, xmid, scal),data=dataToFit,
      # #                      start=list(Asym=startVals[1,1]+200,xmid=startVals[2,1],scal=startVals[3,1]))
      # #  }else if (sid=='8e907') {
      # #    fit.nlme0 <- gnls(PPvec ~ SSlogis(stimLevelVec, Asym, xmid, scal),data=dataToFit,
      # #                      start=list(Asym=startVals[1,1]+200,xmid=startVals[2,1],scal=startVals[3,1]),control = list(maxiter = 500))
      # #  }else {
      # #    fit.nlme0 <- gnls(PPvec ~ SSlogis(stimLevelVec, Asym, xmid, scal),data=dataToFit,
      # #                      start=list(Asym=startVals[1,1]+200,xmid=startVals[2,1],scal=startVals[3,1]))
      # #  }
      # # 
      # #  #Create a range of doses:
      # #  stimLevelVecNew <- data.frame(stimLevelVec = seq(min(dataToFit$stimLevelVec), max(dataToFit$stimLevelVec), length.out = 200))
      # #  #Create a new data frame for ggplot using predict and your range of new
      # #  #doses:
      # #  predData=data.frame(PPvec=predict(fit.nlme0,newdata=stimLevelVecNew),stimLevelVec=stimLevelVecNew)
      # # 
      # #  p3 <- ggplot(dataToFit, aes(x=stimLevelVec, y=PPvec,colour=stimLevelVec)) +
      # #    geom_jitter(width=35) +
      # #    facet_wrap(~blockVec, labeller = as_labeller(blockNames))+
      # #    labs(x = expression(paste("Stimulation current ",mu,"A")),y=expression(paste("Peak to Peak magnitude ",mu,"V"), fill="stimulus level"),title = paste0("Subject ", subjectNum, " ID ", sid," DBS Paired Pulse EP Measurements")) +
      # #    guides(colour=guide_legend("stimulation level"))+ geom_line(data=predData,aes(x=stimLevelVec,y=PPvec),size=2,colour='black')
      # #  print(p3)
      # # 
      # #  if (savePlot && !avgMeas) {
      # #    ggsave(paste0("subj_", subjectNum, "_ID_", sid,"_scatter_sameModel.png"), units="in", width=figWidth, height=figHeight, dpi=600)
      # #  }
      # #  else if (savePlot && avgMeas){
      # #    ggsave(paste0("subj_", subjectNum, "_ID_", sid,"_scatter_sameModel_avg.png"), units="in", width=figWidth, height=figHeight, dpi=600)
      # # 
      # #  }
      # 
      # fit.nlme1 <- nlme(PPvec ~ SSlogis(stimLevelVec, Asym, xmid, scal),
      #                   fixed=list(Asym ~ blockVec, xmid + scal ~ 1),
      #                   random = scal ~ 1|blockVec,
      #                   start=list(fixed=c(Asym=rep(startVals[1,1]+200,length(blockIntLM)),xmid=startVals[2,1],scal=startVals[3,1])),
      #                   data=dataToFit)
      # 
      # # vs.   fixed=list(Asym ~ blockVec, xmid + scal ~ 1),
      # 
      # #Create a range of doses:
      # predDataGroups <-expand.grid(stimLevelVec = seq(min(dataToFit$stimLevelVec), max(dataToFit$stimLevelVec), length.out = 100),blockVec = unique(dataToFit$blockVec))
      # #Create a new data frame for ggplot using predict and your range of new 
      # #doses:
      # predDataGroups$PPvec=predict(fit.nlme1,newdata=predDataGroups)
      # 
      # p4 <- ggplot(dataToFit, aes(x=stimLevelVec, y=PPvec,colour=stimLevelVec)) + 
      #   geom_jitter(width=35) +
      #   facet_wrap(~blockVec, labeller = as_labeller(blockNames))+
      #   #  labs(x = expression(paste("Stimulation current ",mu,"A")),y=expression(paste("Peak to Peak magnitude ",mu,"V"), fill="stimulus level"),title = paste0("Subject ", subjectNum, " ID ", sid," DBS Paired Pulse EP Measurements")) + 
      #   labs(x = expression(paste("Stimulation current ",mu,"A")),y=expression(paste("Peak to Peak magnitude ",mu,"V"), fill="stimulus level"),title = paste0("Baseline and post-conditioning CEP measurements")) + 
      #   guides(colour=guide_colorbar("stimulation level"))+ 
      #   geom_line(data=predDataGroups,aes(x=stimLevelVec,y=PPvec),size=2,colour='black')
      # 
      # print(p4)
      # 
      # if (savePlot && !avgMeas) {
      #   ggsave(paste0("subj_", subjectNum, "_ID_", sid,"_scatter_diffModel.png"), units="in", width=figWidth, height=figHeight, dpi=600)
      # } 
      # else if (savePlot && avgMeas) {
      #   ggsave(paste0("subj_", subjectNum, "_ID_", sid,"_scatter_diffModel_avg.png"), units="in", width=figWidth, height=figHeight, dpi=600)
      # }
      # 
      # xtra ss test
      #   anova(fit.nlme0, fit.nlme1)  
      # Likelihood ratio test
      
      #  lrtest(fit.nlme0, fit.nlme1,fit.glm,fit.lmm)  
      
      # aic
      # AIC(fit.nlme0, fit.nlme1,fit.glm,fit.lmm)  
      
      #BIC(fit.nlme0, fit.nlme1,fit.glm,fit.lmm)
    }
  }
  
  # now do comparison at highest stim level relative to baseline
  dataList <- do.call(rbind,dataList)
  #blockList <- do.call(rbind,blockList)
  
  #plot
  grouped <- group_by(dataList, sidVec, chanVec, blockType,mapStimLevel)
  dataListSummarize <- summarise(grouped,meanPerc = mean(percentDiff),sdPerc = sd(percentDiff),
                                 meanAbs=mean(absDiff), sdDiff=sd(percentDiff),meanPP = mean(PPvec),sdPP = sd(PPvec))
  
  
  p3 <- ggplot(dataList, aes(x=mapStimLevel, y=percentDiff,colour=blockType,fill=blockType)) +
    geom_point(position=position_jitterdodge(dodge.width=0.75)) +geom_smooth(method=lm) +
    labs(x = expression(paste("Stimulation Level")),y=expression(paste("Percent Difference in EP Magnitude from baseline"), fill="stimulus level"),title = paste0("Changes in EP Magnitude"))
  print(p3)
  
  if (savePlot && !avgMeas) {
    ggsave(paste0("across_subj_percent.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    ggsave(paste0("across_subj_percent.eps"), units="in", width=figWidth, height=figHeight, dpi=600,device="eps")
    
  } else if (savePlot && avgMeas){
    ggsave(paste0("across_subj_percent_avg.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    ggsave(paste0("across_subj_percent_avg.eps"), units="in", width=figWidth, height=figHeight, dpi=600,device="eps")
    
  }
  
  p4 <- ggplot(dataList, aes(x=mapStimLevel, y=absDiff,colour=blockType,fill=blockType)) +
    geom_point(position=position_jitterdodge(dodge.width=0.75)) +geom_smooth(method=lm) +
    labs(x = expression(paste("Stimulation Level")),y=expression(paste("Absolute Difference in EP Magnitude from baseline (",mu,"V)"), fill="stimulus level"),title = paste0("Changes in EP Magnitude"))
  print(p4)
  
  if (savePlot && !avgMeas) {
    ggsave(paste0("across_subj_abs.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    ggsave(paste0("across_subj_abs.eps"), units="in", width=figWidth, height=figHeight, dpi=600,device="eps")
    
  } else if (savePlot && avgMeas){
    ggsave(paste0("across_subj_abs_avg.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    ggsave(paste0("across_subj_abs_avg.eps"), units="in", width=figWidth, height=figHeight, dpi=600,device="eps")
    
  }
  
  p5 <- ggplot(dataListSummarize, aes(x=mapStimLevel, y=meanPerc,colour=blockType,fill=blockType)) +
    geom_point(position=position_jitterdodge(dodge.width=0.5)) +geom_smooth(method=lm) +
    labs(x = expression(paste("Stimulation Level")),y=expression(paste("Mean Percent Difference in EP Magnitude"), fill="stimulus level"),title = paste0("Changes in EP Magnitude"))
  print(p5)
  
  if (savePlot && !avgMeas) {
    ggsave(paste0("across_subj_mean_perc.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    ggsave(paste0("across_subj_mean_perc.eps"), units="in", width=figWidth, height=figHeight, dpi=600,device="eps")
    
  } else if (savePlot && avgMeas){
    ggsave(paste0("across_subj_mean_perc_avg.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    ggsave(paste0("across_subj_mean_perc_avg.eps"), units="in", width=figWidth, height=figHeight, dpi=600,device="eps")
    
  }
  
  p6 <- ggplot(dataListSummarize, aes(x=mapStimLevel, y=meanAbs,colour=blockType,fill=blockType)) +
    geom_point(position=position_jitterdodge(dodge.width=0.5)) +geom_smooth(method=lm) +
    labs(x = expression(paste("Stimulation Level")),y=expression(paste("Mean Absolute Difference in EP Magnitude from baseline (",mu,"V)"), fill="stimulus level"),title = paste0("Changes in EP Magnitude"))
  print(p6)
  
  if (savePlot && !avgMeas) {
    ggsave(paste0("across_subj_mean_abs.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    ggsave(paste0("across_subj_mean_abs.eps"), units="in", width=figWidth, height=figHeight, dpi=600,device="eps")
    
  } else if (savePlot && avgMeas){
    ggsave(paste0("across_subj_mean_abs_avg.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    ggsave(paste0("across_subj_mean_abs_avg.eps"), units="in", width=figWidth, height=figHeight, dpi=600,device="eps")
    
  }
  
  p7 <- ggplot(dataList, aes(x=mapStimLevel, y=PPvec,colour=blockType,fill=blockType)) +
    geom_point(position=position_jitterdodge(dodge.width=0.75)) +geom_smooth(method=lm) +
    labs(x = expression(paste("Stimulation Level")),y=expression(paste("Peak to Peak Voltage (",mu,"V)"), fill="stimulus level"),title = paste0("EP Magnitude"))
  print(p7)
  if (savePlot && !avgMeas) {
    ggsave(paste0("across_subj_PP.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    ggsave(paste0("across_subj_PP.eps"), units="in", width=figWidth, height=figHeight, dpi=600,device="eps")
    
  } else if (savePlot && avgMeas){
    ggsave(paste0("across_subj_PP_avg.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    ggsave(paste0("across_subj_PP_avg.eps"), units="in", width=figWidth, height=figHeight, dpi=600,device="eps")
    
  }
  
  
  p8 <- ggplot(dataListSummarize, aes(x=mapStimLevel, y=meanPP,colour=blockType,fill=blockType)) +
    geom_point(position=position_jitterdodge(dodge.width=0.5)) +geom_smooth(method=lm) +
    labs(x = expression(paste("Stimulation Level")),y=expression(paste("Mean Peak to Peak Voltage (",mu,"V)"), fill="stimulus level"),title = paste0("EP Magnitude"))
  print(p8)
  
  if (savePlot && !avgMeas) {
    ggsave(paste0("across_subj_mean_PP.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    ggsave(paste0("across_subj_mean_PP.eps"), units="in", width=figWidth, height=figHeight, dpi=600,device="eps")
    
  } else if (savePlot && avgMeas){
    ggsave(paste0("across_subj_mean_PP_avg.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    ggsave(paste0("across_subj_mean_PP_avg.eps"), units="in", width=figWidth, height=figHeight, dpi=600,device="eps")
    
  }  
  

  
  fit.lmmPP = lmerTest::lmer(PPvec ~ mapStimLevel + blockType + chanInCond + (1|sidVec),data=dataList)

  summary(fit.lmmPP)
  # plot(fit.lm)
  summary(glht(fit.lmmPP,linfct=mcp(blockType="Tukey")))
  summary(glht(fit.lmmPP,linfct=mcp(mapStimLevel="Tukey")))
  
  emmeans(fit.lmmPP, list(pairwise ~ blockType), adjust = "tukey")
  emmeans(fit.lmmPP, list(pairwise ~ mapStimLevel), adjust = "tukey")
  
  emm_s.t <- emmeans(fit.lmmPP, pairwise ~ blockType | mapStimLevel)
  emm_s.t <- emmeans(fit.lmmPP, pairwise ~ mapStimLevel | blockType)
  
  anova(fit.lmmPP)
  tab_model(fit.lmmPP)
  
  fit.lmmdiff = lmerTest::lmer(absDiff ~ mapStimLevel + blockType + chanInCond + (1|sidVec),data=dataList)

  fit.lm = lm(absDiff ~ mapStimLevel + blockType + chanInCond ,data=dataList)
  
  summary(fit.lmmdiff)
  # plot(fit.lm)
  summary(glht(fit.lmmdiff,linfct=mcp(blockType="Tukey")))
  summary(glht(fit.lmmdiff,linfct=mcp(mapStimLevel="Tukey")))
  
  emmeans(fit.lmmdiff, list(pairwise ~ blockType), adjust = "tukey")
  emmeans(fit.lmmdiff, list(pairwise ~ mapStimLevel), adjust = "tukey")
  
  emm_s.t <- emmeans(fit.lmmdiff, pairwise ~ blockType | mapStimLevel)
  emm_s.t <- emmeans(fit.lmmdiff, pairwise ~ mapStimLevel | blockType)
  
  anova(fit.lmmdiff)
  tab_model(fit.lmmdiff)
  
  
  
  
}

