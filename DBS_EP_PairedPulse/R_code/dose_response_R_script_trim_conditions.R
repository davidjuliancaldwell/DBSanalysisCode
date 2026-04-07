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
library('effectsize')
library('caret')
library('DescTools')
#library('glmmTMB') # no longer used after switching to lmer

rootDir = here()
dataDir = here("DBS_EP_PairedPulse","R_data")
codeDir = here("DBS_EP_PairedPulse","R_code")
outputDir = here("DBS_EP_PairedPulse","R_output")
dir.create(outputDir, showWarnings = FALSE)

sidVec <- c('46c2a','c963f','2e114','fe7df','e6f3c','9f852',
            '8e907','08b13','e9c9b','41a73','68574',
            '01fee','a23ed')

# NOTE: 46c2a classified as MD (verify against clinical records)
diseaseVec <- c('MD','PD','MD','PD','PD','PD','PD','MD',
                'PD','PD','PD','PD','MD')
repeatedMeasures = TRUE # if true, does repeated measures analysis, if false, does more of ANCOVA style analysis
log_data = TRUE
box_data = FALSE
trim_data = FALSE
min_stim_level = 1
savePlot = 1
avgMeasVec = c(0)
figWidth = 8 
figHeight = 6 

oldwd <- getwd()
setwd(outputDir)

for (avgMeas in avgMeasVec) {

  dataList = list()
  blockList = list()
  conditionList = list()
  index = 1
  
  subjectNumIndex = 1
  
  for (sid in sidVec){
    source(here("DBS_EP_PairedPulse","R_config_files",paste0("subj_",sid,".R")))
    
    brodmann_areas <- read.csv(here("DBS_EP_PairedPulse","R_config_files",paste0(sid,'_MNIcoords_labelled.csv')),header=TRUE,sep = ",",stringsAsFactors=F)
    dataPP <- read.table(here("DBS_EP_PairedPulse","R_data",paste0(sid,'_PairedPulseData_avg_5.csv')),header=TRUE,sep = ",",stringsAsFactors=F, colClasses=c("stimLevelVec"="numeric","sidVec"="factor"))

    if (sid == "2e114"){
      goodVec = c(1000,2500,3500)
      dataPP <- dataPP %>% filter(stimLevelVec %in% goodVec)
    }

    # PPvec is in volts, convert to microvolts
    dataPP$PPvec = dataPP$PPvec*1e6
    dataPP$stimLevelVec = dataPP$stimLevelVec/1e3
    

    # denote which channel was in the conditioning pair 
    
    dataPP <- subset(dataPP,(chanVec %in% chanIntVec))
    
    dataPP$chanInCond <- mapvalues(dataPP$chanVec,
                                   from=chanIntVec,
                                   to=chanIntConditioningPairVec)
    
    dataPP$chanInCond <- as.factor(dataPP$chanInCond)
    # change to factor 
    dataPP$blockVec = as.factor(dataPP$blockVec)
    dataPP$chanVec = as.factor(dataPP$chanVec)
    print(sid)
    
    dataPP$disease = as.factor(diseaseVec[subjectNumIndex])
    dataPP$subjectNum = as.factor(subjectNumIndex)
    
    
    for (chanInt in chanIntVec){
      
      # select data of interest
      dataInt <- subset(dataPP,(chanVec == chanInt) & (blockVec %in% blockIntPlot))
      dataInt <- na.exclude(dataInt)
      dataInt$blockVec <- droplevels(dataInt$blockVec)
      
      # map linear subject numbering, rather than using the total subject # order from each subj_ setup file
      
      
      # add disease
      
      # map stimulation levels to consistent ordering for between subject comparison
      uniqueStimLevel = as.double(unique(dataInt$stimLevelVec))
      
      if (sid == "2e114"){
        mappingStimLevel =c(1,3,4)
        
      } else if (sid == "41a73" | sid == "68574"){
        mappingStimLevel = c(1,3,4,0)
      }
      else {
        mappingStimLevel =c(1:length(uniqueStimLevel))
      }
      dataInt$mapStimLevel<- mapvalues(dataInt$stimLevelVec,
                                        from=sort(uniqueStimLevel),
                                        to=mappingStimLevel)
      
     # uniqueBlockLevel = unique(dataInt$blockVec)
     # blockTypeTrim = blockType[uniqueBlockLevel]
      
     # dataInt$blockType <- mapvalues(dataInt$blockVec,
      #                               from=uniqueBlockLevel,
       #                              to=blockTypeTrim)
      
      
      dataInt$blockType <- mapvalues(dataInt$blockVec,
                                     from=blockIntPlot,
                                    to=blockType)
      
      dataInt$mapStimLevel = as.ordered(dataInt$mapStimLevel)
      dataInt$blockType = as.factor(dataInt$blockType)
      
      # get which brodmann area label is for this electrode
      dataInt$baLabel = as.character(brodmann_areas$ba.label[chanInt])
      dataInt$aalLabel = as.character(brodmann_areas$aal.label[chanInt])
      
      # get all to same side
      dataInt$baLabel[dataInt$baLabel=='Right-PrimMotor (4)'] = 'Left-PrimMotor (4)'
      dataInt$aalLabel[dataInt$aalLabel=='Postcentral_R'] = 'Postcentral_L'
      
      dataInt$baLabel[is.na(dataInt$baLabel)] <- "Unknown"
      dataInt$aalLabel[is.na(dataInt$aalLabel)] <- "Unknown"
      dataInt$baLabel = as.factor(dataInt$baLabel)
      dataInt$aalLabel = as.factor(dataInt$aalLabel)
 
      dataInt <- subset(dataInt, PPvec<1000)
      dataInt <- subset(dataInt, PPvec>10)
      
      dataInt <- na.exclude(dataInt)
      
      if (trim_data){
        dataInt <- dataInt %>% group_by(blockVec,blockType,mapStimLevel) %>% mutate(PPvecLabel = !is.element(seq_len(length(PPvec)),attr(Trim(PPvec,0.025,na.rm=FALSE),'trim')))
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
      

      # Assign first/second half based on trial order within each (block, stim level).
      # Row order in the CSV preserves temporal sequence from MATLAB.
      dataInt <- dataInt %>%
        group_by(blockVec, mapStimLevel) %>%
        mutate(halfBlock = ifelse(row_number() <= ceiling(n()/2), "first", "second")) %>%
        ungroup()
      dataInt$halfBlock <- as.factor(dataInt$halfBlock)

      for (comparison in whichCompareVec){
        dataIntCompare <- subset(dataInt,(blockVec %in% comparison))

        #if(!(is.data.frame(dataIntCompare) & (nrow(dataIntCompare)==0))){

        # Filter per stim level: require >= 3 trials per (stim level, block)
        # cell, and both blocks present at that stim level. Stim levels that
        # fail are dropped individually rather than excluding the entire subject.
        dataIntCompare <- dataIntCompare %>%
          group_by(mapStimLevel, blockVec) %>% mutate(enough = n() >= 3) %>%
          group_by(mapStimLevel) %>% filter(all(enough) & length(unique(blockVec)) > 1) %>%
          ungroup()
        
        if(!(empty(dataIntCompare))){   
          
        dataIntCompare <- dataIntCompare %>%
          group_by(mapStimLevel) %>%
          mutate(
                  absDiff = PPvec - mean(PPvec[blockVec==comparison[1]]),
                 percentDiff = 100*(PPvec - mean(PPvec[blockVec==comparison[1]]))/mean(PPvec[blockVec==comparison[1]]) )
        
        dataIntCompare <- dataIntCompare %>% 
          group_by(mapStimLevel,chanVec,blockType) %>% 
          mutate(
            meanPP = mean(PPvec) ) 
        
        dataIntCompare = as_data_frame(dataIntCompare)
        dataIntCompare$index = index
        
        name <- unique(dataIntCompare %>% filter(blockType!='baseline')%>%select(blockType))
        name <- as.character(name$blockType)
        
        dataIntCompare$overallBlockType <- name
        dataIntCompare <- dataIntCompare %>% mutate(pre_post = case_when(blockType == name ~ 'post',
          blockType=='baseline' ~'pre'
        ))
        
        dataList[[index]] = dataIntCompare
        blockList[[index]] = comparison
        
        #fit.lm    = lm(PPvec ~ mapStimLevel + blockVec + mapStimLevel*blockVec,data=dataIntCompare)
        #fit.lm    = lm(PPvec ~ mapStimLevel + blockVec,data=dataIntCompare)
        
        #summary(fit.lm)
        # plot(fit.lm)
        #summary(glht(fit.lm,linfct=mcp(blockVec="Tukey")))
        #summary(glht(fit.lm,linfct=mcp(mapStimLevel="Tukey")))
        
        #emmeans(fit.lm, list(pairwise ~ blockVec), adjust = "tukey")
        #emmeans(fit.lm, list(pairwise ~ mapStimLevel), adjust = "tukey")
        
        #emm_s.t <- emmeans(fit.lm, pairwise ~ blockVec | mapStimLevel)
        #emm_s.t <- emmeans(fit.lm, pairwise ~ mapStimLevel | blockVec)
        
        #anova(fit.lm)
        #tab_model(fit.lm)
        
        index = index + 1
        }
        
      }
      
      subjectNumInterest = unique(dataInt$subjectNum)
      
      tryCatch({
        p <- ggplot(dataInt, aes(x=stimLevelVec, y=PPvec,color=stimLevelVec)) +
          geom_point(position=position_jitter(width=0.1)) +  geom_smooth(method=lm) + facet_wrap(~blockVec, labeller = as_labeller(blockNames))+
          labs(x = expression(paste("Stimulation Current (mA)")),y=expression(paste("Peak to Peak Magnitude (",mu,"V)")),title = paste0("Subject ", subjectNumInterest," DBS Paired Pulse EP Measurements")) +
          guides(colour=guide_colorbar("Stimulation Level"))
        print(p)

        if (savePlot && !avgMeas) {
          ggsave(paste0("subj_", subjectNum, "_ID_", sid,"_chan_",chanInt,"_scatter_lm_AUC.png"), units="in", width=figWidth, height=figHeight, dpi=600)
          ggsave(paste0("subj_", subjectNum, "_ID_", sid,"_chan_",chanInt,"_scatter_lm_AUC.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
        }
        else if (savePlot && avgMeas){
          ggsave(paste0("subj_", subjectNum, "_ID_", sid,"_chan_",chanInt,"_scatter_lm_avg_AUC.png"), units="in", width=figWidth, height=figHeight, dpi=600)
          ggsave(paste0("subj_", subjectNum, "_ID_", sid,"_chan_",chanInt,"_scatter_lm_avg_AUC.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
        }
      }, error=function(e) { cat("  Skipping per-block scatter for", sid, "chan", chanInt, ":", e$message, "\n") })
      
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
      # #    ggsave(paste0("subj_", subjectNum, "_ID_", sid,"_scatter_sameModel_AUC.png"), units="in", width=figWidth, height=figHeight, dpi=600)
      # #  }
      # #  else if (savePlot && avgMeas){
      # #    ggsave(paste0("subj_", subjectNum, "_ID_", sid,"_scatter_sameModel_avg_AUC.png"), units="in", width=figWidth, height=figHeight, dpi=600)
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
      #   ggsave(paste0("subj_", subjectNum, "_ID_", sid,"_scatter_diffModel_AUC.png"), units="in", width=figWidth, height=figHeight, dpi=600)
      # } 
      # else if (savePlot && avgMeas) {
      #   ggsave(paste0("subj_", subjectNum, "_ID_", sid,"_scatter_diffModel_avg_AUC.png"), units="in", width=figWidth, height=figHeight, dpi=600)
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
    subjectNumIndex = subjectNumIndex + 1
  }
  
  # now do comparison at highest stim level relative to baseline
  dataList <- do.call(rbind,dataList)
  #blockList <- do.call(rbind,blockList)
  
  if (repeatedMeasures){
    # Include A/A 25 (same-polarity 25ms control, 3 subjects) alongside
    # A/B 25, A/B 200, and A/A 200. A/B 100 excluded (1 subject only).
    dataList <- dataList %>% filter((blockType %in% c('baseline','A/B 25','A/B 200','A/A 200','A/A 25')) & (overallBlockType %in% c('A/B 25','A/B 200','A/A 200','A/A 25')))
    dataList$overallBlockType <- as.factor(dataList$overallBlockType)
  }
  else {
    dataList <- dataList %>% filter(blockType %in% c('baseline','A/B 25','A/B 200','A/A 200','A/A 25'))
  }

  cat("\n=== Subjects in analysis ===\n")
  cat("Unique subjects:", length(unique(dataList$subjectNum)), "\n")
  cat("Subject IDs:", paste(unique(dataList$sidVec), collapse=", "), "\n")
  cat("Conditions:", paste(unique(dataList$overallBlockType), collapse=", "), "\n")
  cat("Per-subject obs:", paste(tapply(dataList$PPvec, dataList$sidVec, length), collapse=", "), "\n")

  # Aggregate to cell medians to avoid trial-level pseudoreplication.
  # Each row becomes the median across trials within a
  # (subject, channel, block, stim-level) cell. Median is robust to
  # outliers without requiring trimming (trim_data = FALSE).
  # Cell medians without halfBlock split (for primary model)
  dataListAgg <- dataList %>%
    group_by(subjectNum, chanVec, mapStimLevel, blockVec, blockType,
             overallBlockType, pre_post, disease, chanInCond, baLabel,
             aalLabel, sidVec, index) %>%
    summarise(PPvec = median(PPvec),
              absDiff = median(absDiff),
              percentDiff = median(percentDiff),
              meanPP = first(meanPP),
              stimLevelVec = first(stimLevelVec),
              n_trials = n(),
              .groups = "drop")

  # Cell medians WITH halfBlock split (for secondary halfBlock model)
  dataListAggHalf <- dataList %>%
    group_by(subjectNum, chanVec, mapStimLevel, blockVec, blockType,
             overallBlockType, pre_post, disease, chanInCond, baLabel,
             aalLabel, sidVec, index, halfBlock) %>%
    summarise(PPvec = median(PPvec),
              absDiff = median(absDiff),
              percentDiff = median(percentDiff),
              meanPP = first(meanPP),
              stimLevelVec = first(stimLevelVec),
              n_trials = n(),
              .groups = "drop")

  #plot
  grouped <- group_by(dataList, sidVec, chanVec, blockType,mapStimLevel,disease)
  dataListSummarize <- summarise(grouped,meanPerc = mean(percentDiff),sdPerc = sd(percentDiff),
                                 meanAbs=mean(absDiff), sdDiff=sd(percentDiff),meanPP = mean(PPvec),sdPP = sd(PPvec))
  
  p3 <- ggplot(dataList, aes(x=as.numeric(mapStimLevel), y=percentDiff,color=blockType)) +
    geom_point(position=position_jitterdodge(dodge.width=0.75)) +geom_smooth(method=lm) +
    labs(x = expression(paste("Stimulation Level")),y=expression(paste("Percent Difference in EP Magnitude from Baseline")), color="Experimental Condition",title = paste0("Changes in EP Magnitude"))
  print(p3)
  
  if (savePlot && !avgMeas) {
    ggsave(paste0("across_subj_percent_AUC.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    ggsave(paste0("across_subj_percent_AUC.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
    
  } else if (savePlot && avgMeas){
    ggsave(paste0("across_subj_percent_avg_AUC.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    ggsave(paste0("across_subj_percent_avg_AUC.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
    
  }
  
  p4 <- ggplot(dataList, aes(x=as.numeric(mapStimLevel), y=absDiff,color=blockType)) +
    geom_point(position=position_jitterdodge(dodge.width=0.75)) +geom_smooth(method=lm) +
    labs(x = expression(paste("Stimulation Level")),y=expression(paste("Absolute Difference in EP Magnitude from Baseline (",mu,"V)")),color="Experimental Condition",title = paste0("Changes in EP Magnitude"))
  print(p4)
  
  if (savePlot && !avgMeas) {
    ggsave(paste0("across_subj_abs_AUC.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    ggsave(paste0("across_subj_abs_AUC.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
    
  } else if (savePlot && avgMeas){
    ggsave(paste0("across_subj_abs_avg_AUC.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    ggsave(paste0("across_subj_abs_avg_AUC.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
    
  }
  
  p5 <- ggplot(dataListSummarize, aes(x=as.numeric(mapStimLevel), y=meanPerc,color=blockType)) +
    geom_point(position=position_jitterdodge(dodge.width=0.5)) +geom_smooth(method=lm) +
    labs(x = expression(paste("Stimulation Level")),y=expression(paste("Mean Percent Difference in EP Magnitude")), color="Experimental Condition",title = paste0("Changes in EP Magnitude"))
  print(p5)
  
  if (savePlot && !avgMeas) {
    ggsave(paste0("across_subj_mean_perc_AUC.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    ggsave(paste0("across_subj_mean_perc_AUC.eps"), units="in", width=figWidth, height=figHeight, device=cairo_ps, fallback_resolution=600)
    
  } else if (savePlot && avgMeas){
    ggsave(paste0("across_subj_mean_perc_avg_AUC.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    ggsave(paste0("across_subj_mean_perc_avg_AUC.eps"), units="in", width=figWidth, height=figHeight, device=cairo_ps, fallback_resolution=600)
    
  }
  
  p6 <- ggplot(dataListSummarize, aes(x=mapStimLevel, y=meanAbs,color=blockType)) +
    geom_point(position=position_jitterdodge(dodge.width=0.5)) +geom_smooth(method=lm) +
    labs(x = expression(paste("Stimulation Level")),y=expression(paste("Mean Absolute Difference in EP Magnitude from Baseline (",mu,"V)")),
         color="Experimental Condition",title = paste0("Changes in EP Magnitude"))
  print(p6)
  
  if (savePlot && !avgMeas) {
    ggsave(paste0("across_subj_mean_abs_AUC.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    ggsave(paste0("across_subj_mean_abs_AUC.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
    
  } else if (savePlot && avgMeas){
    ggsave(paste0("across_subj_mean_abs_avg_AUC.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    ggsave(paste0("across_subj_mean_abs_avg_AUC.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
    
  }
  
  p7 <- ggplot(dataList, aes(x=mapStimLevel, y=PPvec,color=blockType)) +
    geom_point(position=position_jitterdodge(dodge.width=0.75)) +geom_smooth(method=lm) +
    labs(x = expression(paste("Stimulation Level")),y=expression(paste("Peak to Peak Voltage (",mu,"V)")),
         color="Experimental Condition",title = paste0("EP Magnitude"))
  print(p7)
  if (savePlot && !avgMeas) {
    ggsave(paste0("across_subj_PP_AUC.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    ggsave(paste0("across_subj_PP_AUC.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
    
  } else if (savePlot && avgMeas){
    ggsave(paste0("across_subj_PP_avg_AUC.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    ggsave(paste0("across_subj_PP_avg_AUC.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
    
  }

  p8 <- ggplot(dataListSummarize, aes(x=mapStimLevel, y=meanPP,color=blockType)) +
    geom_point(position=position_jitterdodge(dodge.width=0.5)) +geom_smooth(method=lm) +
    labs(x = expression(paste("Stimulation Level")),y=expression(paste("Mean peak to peak voltage (",mu,"V)")), color="Experimental condition",title = paste0("EP magnitude"))
  print(p8)
  
  if (savePlot && !avgMeas) {
    ggsave(paste0("across_subj_mean_PP_AUC.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    ggsave(paste0("across_subj_mean_PP_AUC.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
    
  } else if (savePlot && avgMeas){
    ggsave(paste0("across_subj_mean_PP_avg_AUC.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    ggsave(paste0("across_subj_mean_PP_avg_AUC.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
    
  }  
  
  # figWidth = 8
  # figHeight = 4
  # p11 <- ggplot(data = dataListSummarize, aes(x = blockType, y = meanAbs,color=blockType)) +
  #   geom_boxplot(notch=TRUE,outlier.shape=NA)  + geom_jitter(shape=16, position=position_jitter(0.2),aes(alpha = mapStimLevel)) +
  #   labs(x = expression(paste("Experimental Condition")),y=expression(paste("Absolute Difference from Baseline Peak-To-Peak (",mu,"V)")),color="Experimental Condition",alpha="Ordered Stim Level",title = paste0("EP Difference from Baseline by Conditioning Protocol"))
  # print(p11)
  # 
  # if (savePlot && !avgMeas) {
  #   ggsave(paste0("across_subj_mean_abs_box_AUC.png"), units="in", width=figWidth, height=figHeight, dpi=600)
  #   ggsave(paste0("across_subj_mean_abs_box_AUC.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
  # 
  # } else if (savePlot && avgMeas){
  #   ggsave(paste0("across_subj_mean_abs_box_avg_AUC.png"), units="in", width=figWidth, height=figHeight, dpi=600)
  #   ggsave(paste0("across_subj_mean_abs_box_avg_AUC.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
  # 
  # }
  
  # box plot
  
  figWidth = 8
  figHeight = 4
  goodVecBlock <- c("A/B 25","A/B 200","A/A 200","A/A 25")
  p11 <- ggplot(data =  dataListSummarize%>% filter(blockType %in% goodVecBlock), aes(x = blockType, y = meanAbs,color=blockType)) +
    geom_boxplot(notch=TRUE,outlier.shape=NA)  + geom_jitter(shape=16, position=position_jitter(0.2),aes(alpha = mapStimLevel)) +
    labs(x = expression(paste("Experimental Condition")),y=expression(paste("Absolute Difference from Baseline Peak-To-Peak (",mu,"V)")),color="Experimental Condition",alpha="Ordered Stim Level",title = paste0("EP Difference from Baseline by Conditioning Protocol")) +
    scale_color_brewer(palette="Dark2")
  print(p11)
  
  if (savePlot && !avgMeas) {
    ggsave(paste0("across_subj_mean_abs_box_all_subjs_AUC.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    ggsave(paste0("across_subj_mean_abs_box_all_subjs_AUC.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
    
  } else if (savePlot && avgMeas){
    ggsave(paste0("across_subj_mean_abs_box_avg_AUC.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    ggsave(paste0("across_subj_mean_abs_box_avg_AUC.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
    
  }
  
  # box plot — match color order and dimensions with effect size plot

  goodVecBlock <- c("A/A 200","A/A 25","A/B 200","A/B 25")
  pct_data <- dataListSummarize %>% filter(blockType %in% goodVecBlock) %>%
    mutate(blockType = factor(blockType, levels=goodVecBlock))
  p11 <- ggplot(data = pct_data, aes(x = blockType, y = meanPerc, color=blockType)) +
    geom_boxplot(notch=TRUE,outlier.shape=NA) +
    geom_jitter(shape=16, position=position_jitter(0.2), size=2.5, aes(alpha = mapStimLevel)) +
    labs(x = "Experimental Condition",
         y = "Percent Difference from Matched Baseline (log scale)",
         color="Condition", alpha="Ordered Stim Level",
         title = "Percent Change from Baseline by Conditioning Protocol") +
    scale_color_brewer(palette="Dark2") +
    theme_bw(base_size=16) +
    theme(plot.title = element_text(size=16, face="bold"),
          axis.title = element_text(size=14),
          axis.text = element_text(size=12),
          legend.title = element_text(size=12),
          legend.text = element_text(size=11))
  print(p11)

  if (savePlot && !avgMeas) {
    ggsave(paste0("across_subj_pct_diff_baseline_box_all_subjs_AUC.png"), units="in", width=7, height=7, dpi=600)
    ggsave(paste0("across_subj_pct_diff_baseline_box_all_subjs_AUC.eps"), units="in", width=7, height=7, device=cairo_ps, fallback_resolution=600)
    
  } else if (savePlot && avgMeas){
    ggsave(paste0("across_subj_pct_diff_baseline_box_avg_AUC.png"), units="in", width=7, height=7, dpi=600)
    ggsave(paste0("across_subj_pct_diff_baseline_box_avg_AUC.eps"), units="in", width=7, height=7, device=cairo_ps, fallback_resolution=600)
  }

  # --- Plot: Cell-median percent difference from matched baseline ---
  # Compute percent difference directly from cell-median PPvec values.
  # Each post cell is matched to its baseline cell by (subject, channel,
  # stim level, comparison index). index uniquely identifies each
  # baseline-conditioning pair, so the join gives exactly one baseline
  # per post cell.
  baseline_medians <- dataListAgg %>%
    filter(pre_post == "pre") %>%
    dplyr::select(subjectNum, chanVec, mapStimLevel, index, PPvec_baseline = PPvec)
  post_medians <- dataListAgg %>%
    filter(pre_post == "post") %>%
    left_join(baseline_medians, by = c("subjectNum","chanVec","mapStimLevel","index"))
  post_medians$pctDiffMedian <- 100 * (post_medians$PPvec - post_medians$PPvec_baseline) / post_medians$PPvec_baseline

  goodVecBlock_model <- c("A/A 200","A/A 25","A/B 200","A/B 25")
  model_pct_data <- post_medians %>%
    filter(overallBlockType %in% goodVecBlock_model) %>%
    mutate(overallBlockType = factor(overallBlockType, levels=goodVecBlock_model))
  p_model_data <- ggplot(data = model_pct_data, aes(x = overallBlockType, y = pctDiffMedian, color = overallBlockType)) +
    geom_boxplot(notch=TRUE, outlier.shape=NA) +
    geom_jitter(shape=16, position=position_jitter(0.2), size=2.5, aes(alpha = mapStimLevel)) +
    geom_hline(yintercept=0, linetype="dashed", color="grey50") +
    labs(x = "Experimental Condition",
         y = "Percent Difference from Matched Baseline (log scale)",
         color="Condition", alpha="Ordered Stim Level",
         title = "Percent Change from Baseline (Cell Medians)") +
    scale_color_brewer(palette="Dark2") +
    theme_bw(base_size=16) +
    theme(plot.title = element_text(size=16, face="bold"),
          axis.title = element_text(size=14),
          axis.text = element_text(size=12),
          legend.title = element_text(size=12),
          legend.text = element_text(size=11))
  print(p_model_data)

  emmFigSize = 7
  if (savePlot && !avgMeas) {
    ggsave(plot=p_model_data, paste0("across_subj_median_log_difference_model_data.png"), units="in", width=emmFigSize, height=emmFigSize, dpi=600)
    ggsave(plot=p_model_data, paste0("across_subj_median_log_difference_model_data.eps"), units="in", width=emmFigSize, height=emmFigSize, device=cairo_ps, fallback_resolution=600)
  }
  
  dataBarInt = dataList%>% filter(blockType %in% goodVecBlock)
  groupedAll <- group_by(dataBarInt,blockType)
  dataListSummarizeBar <- summarise(groupedAll,meanPerc = mean(percentDiff),sdPerc = sd(percentDiff),
                                    meanAbs=mean(absDiff),sdAbs=sd(absDiff), sdDiff=sd(percentDiff),meanPP = mean(PPvec),sdPP = sd(PPvec),
                                    n=n(),seAbs=sdAbs/sqrt(n))
  
  # bar plot with error
  figWidth = 8
  figHeight = 4
  goodVecBlock <- c("A/B 25","A/B 200","A/A 200","A/A 25")
  p12 <- ggplot(data =  dataListSummarizeBar, aes(x = blockType, y = meanAbs,fill=blockType)) +
    geom_bar(stat="identity")  +
    labs(x = expression(paste("Experimental Condition")),y=expression(paste("Absolute Difference from Baseline Peak-To-Peak (",mu,"V)")),fill="Experimental Condition",title = paste0("EP Difference from Baseline by Conditioning Protocol")) +
    scale_fill_brewer(palette="Dark2") +
    geom_errorbar(aes(ymin=meanAbs-seAbs, ymax=meanAbs+seAbs), width=.2,
                  position=position_dodge(.9)) 
  print(p12)
  
  if (savePlot && !avgMeas) {
    ggsave(paste0("across_subj_mean_abs_bar_standarderror_AUC.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    ggsave(paste0("across_subj_mean_abs_bar_standarderror_AUC.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
    
  } else if (savePlot && avgMeas){
    ggsave(paste0("across_subj_mean_abs_bar_standarderror_avg_AUC.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    ggsave(paste0("across_subj_mean_abs_bar_standarderror_avg_AUC.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
    
  }
  
  
  # bar plot 
  figWidth = 8
  figHeight = 4
  goodVecBlock <- c("A/B 25","A/B 200","A/A 200","A/A 25")
  p14 <- ggplot(data =  dataListSummarizeBar, aes(x = blockType, y = meanAbs,fill=blockType)) +
    geom_bar(stat="identity")  +
    labs(x = expression(paste("Experimental Condition")),y=expression(paste("Absolute Difference from Baseline Peak-To-Peak (",mu,"V)")),fill="Experimental Condition",title = paste0("EP Difference from Baseline by Conditioning Protocol")) +
    scale_fill_brewer(palette="Dark2")
  print(p14)
  
  if (savePlot && !avgMeas) {
    ggsave(paste0("across_subj_mean_abs_bar_AUC.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    ggsave(paste0("across_subj_mean_abs_bar_AUC.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
    
  } else if (savePlot && avgMeas){
    ggsave(paste0("across_subj_mean_abs_bar_avg_AUC.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    ggsave(paste0("across_subj_mean_abs_bar_avg_AUC.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
    
  }
  
  
  # bar with dots
  figWidth = 8
  figHeight = 4
  goodVecBlock <- c("A/B 25","A/B 200","A/A 200","A/A 25")
  p13 <- ggplot() +
    geom_jitter(data =  dataListSummarize%>% filter(blockType %in% goodVecBlock),shape=16, position=position_jitter(0.2), aes(x = blockType, y = meanAbs,color=blockType,alpha = mapStimLevel)) +
    guides(fill=FALSE) +
    geom_bar(data=dataListSummarizeBar,stat="identity",aes(x = blockType, y = meanAbs,fill=blockType))  +
    labs(x = expression(paste("Experimental Condition")),y=expression(paste("Absolute Difference from Baseline Peak-To-Peak (",mu,"V)")),color="Experimental Condition",alpha="Ordered Stim Level",title = paste0("EP Difference from Baseline by Conditioning Protocol")) +
    scale_color_brewer(palette="Dark2") +
    scale_fill_brewer(palette="Dark2")
  print(p13)
  
  if (savePlot && !avgMeas) {
    ggsave(paste0("across_subj_mean_abs_bar_dots_AUC.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    ggsave(paste0("across_subj_mean_abs_bar_dots_AUC.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
    
  } else if (savePlot && avgMeas){
    ggsave(paste0("across_subj_mean_abs_bar_dots_avg_AUC.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    ggsave(paste0("across_subj_mean_abs_bar_dots_avg_AUC.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
    
  }

  
  # all subjects
  figHeight = 10
  figWidth = 8
  p9 <- ggplot(dataList, aes(x=as.numeric(stimLevelVec), y=absDiff,color=blockType)) +facet_wrap(~subjectNum,scales = "free") +
    geom_point(position=position_jitterdodge(dodge.width=0.5)) +geom_smooth(method=lm) +
    labs(x = expression(paste("Stimulation Level (mA)")),y=expression(paste("Absolute Difference in EP Magnitude from Baseline (",mu,"V)")),color="Experimental Condition",title = paste0("Changes in EP magnitude")) +
    theme(legend.position = c(0.9, 0.1))
  print(p9)
  
  if (savePlot && !avgMeas) {
    ggsave(paste0("diff_scale_each_subj_abs_AUC.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    ggsave(paste0("diff_scale_each_subj_abs_AUC.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
    
  } else if (savePlot && avgMeas){
    ggsave(paste0("diff_scale_each_subj_abs_avg_AUC.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    ggsave(paste0("diff_scale_each_subj_abs_avg_AUC.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
    
  }
  if (repeatedMeasures){

  # === PRIMARY MODEL: Cell-mean aggregation with crossed RE + channel slope ===
  # Cell means eliminate within-cell pseudoreplication (non-independent
  # trials per cell collapsed to 1 value).
  # Crossed RE structure because multi-channel subjects (6, 7, 11, 12) had
  # channels recorded simultaneously — blocks are crossed with channels,
  # not nested within them.
  #   (1|subjectNum)              — between-subject variance (9 groups)
  #   (1+stim_linpoly|subjectNum:chanVec) — channel-level intercept + linear
  #       dose-response slope (13 groups). Random slope allows each electrode
  #       to have its own dose-response steepness, reflecting cortical
  #       position. Model comparison showed channel slope (AIC=158) beats
  #       block slope (AIC=210) and both-slopes (AIC=161, block slope
  #       collapses to ~0). Dose-response steepness is a property of
  #       electrode location, not recording session.
  #   (1|subjectNum:blockVec)     — block-level intercept (41 groups)
  # Disease and baLabel dropped (underpowered: 9 subjects, 13 channels).
  #
  # stim_linpoly: linear polynomial contrast extracted from contr.poly().
  # This is the .L component of the ordered mapStimLevel factor, provided
  # as a numeric column so it can be used as a random slope.
  n_stim_levels <- nlevels(dataListAgg$mapStimLevel)
  poly_lin <- contr.poly(n_stim_levels)[, 1]
  dataListAgg$stim_linpoly <- poly_lin[as.numeric(dataListAgg$mapStimLevel)]
  dataList$stim_linpoly <- poly_lin[as.numeric(dataList$mapStimLevel)]

  fit.lmmPP = lmerTest::lmer(PPvec ~ mapStimLevel + chanInCond +
    overallBlockType*pre_post +
    (1 + stim_linpoly|subjectNum:chanVec) + (1|subjectNum) + (1|subjectNum:blockVec),
    data=dataListAgg)

  # --- halfBlock model: three-way interaction ---
  # Tests whether the conditioning effect (overallBlockType x pre_post)
  # differs between the first and second half of trials within each block.
  # Uses dataListAggHalf (split by halfBlock) rather than dataListAgg.
  dataListAggHalf$stim_linpoly <- poly_lin[as.numeric(dataListAggHalf$mapStimLevel)]
  cat("\n=== halfBlock model: three-way interaction ===\n")
  cat("Observations (split):", nrow(dataListAggHalf), "vs primary:", nrow(dataListAgg), "\n")
  fit.lmmPP.half = lmerTest::lmer(PPvec ~ mapStimLevel + chanInCond +
    overallBlockType*pre_post*halfBlock +
    (1 + stim_linpoly|subjectNum:chanVec) + (1|subjectNum) + (1|subjectNum:blockVec),
    data=dataListAggHalf)
  cat("Singular:", isSingular(fit.lmmPP.half), "\n")
  cat("Observations:", nobs(fit.lmmPP.half), "\n")
  print(summary(fit.lmmPP.half))
  cat("\n=== ANOVA (halfBlock model) ===\n")
  print(anova(fit.lmmPP.half, type=3))
  cat("\n=== emmeans: pre_post x halfBlock | overallBlockType ===\n")
  emm_half <- emmeans(fit.lmmPP.half, ~ pre_post * halfBlock | overallBlockType)
  print(contrast(emm_half, interaction="pairwise"))

  # --- Comparison model: Cell-mean, block slope (not channel slope) ---
  fit.lmmPP.blockslope = lmerTest::lmer(PPvec ~ mapStimLevel + chanInCond +
    overallBlockType*pre_post +
    (1 + stim_linpoly|subjectNum:blockVec) + (1|subjectNum/chanVec),
    data=dataListAgg)

  # --- Comparison model: Cell-mean, intercept-only (no random slope) ---
  fit.lmmPP.intonly = lmerTest::lmer(PPvec ~ mapStimLevel + chanInCond +
    overallBlockType*pre_post + (1|subjectNum/chanVec) + (1|subjectNum:blockVec),
    data=dataListAgg)

  # --- Comparison model: Trial-level, flat block RE ---
  fit.lmmPP.trial.flat = lmerTest::lmer(PPvec ~ mapStimLevel + chanInCond +
    overallBlockType*pre_post + (1|subjectNum:chanVec:blockVec), data=dataList)

  # --- Compare approaches ---
  cat("\n=== PRIMARY: Cell-mean, crossed RE + channel slope ===\n")
  cat("Observations:", nobs(fit.lmmPP), "\n")
  print(summary(fit.lmmPP))
  cat("\n=== Comparison: Cell-mean, block slope ===\n")
  cat("Observations:", nobs(fit.lmmPP.blockslope), "\n")
  print(summary(fit.lmmPP.blockslope))
  cat("\n=== Comparison: Cell-mean, intercept-only ===\n")
  cat("Observations:", nobs(fit.lmmPP.intonly), "\n")
  print(summary(fit.lmmPP.intonly))
  cat("\n=== Comparison: Trial-level, flat block RE ===\n")
  cat("Observations:", nobs(fit.lmmPP.trial.flat), "\n")
  print(summary(fit.lmmPP.trial.flat))
  cat("\n=== AIC comparison ===\n")
  cat("Channel slope (primary): ", AIC(fit.lmmPP), "\n")
  cat("Block slope:             ", AIC(fit.lmmPP.blockslope), "\n")
  cat("Intercept-only:          ", AIC(fit.lmmPP.intonly), "\n")
  cat("Trial-level flat RE:     ", AIC(fit.lmmPP.trial.flat), "\n")

  # Use Satterthwaite df (fast on 84 obs). For t_crit used in CI plots
  # and effect sizes, use the median Satterthwaite df across block-level
  # fixed effects as a representative value.
  emm_options(lmer.df = "satterthwaite")
  block_level_df <- summary(fit.lmmPP)$coefficients[
    c("overallBlockTypeA/B 200","overallBlockTypeA/B 25",
      "pre_postpre","overallBlockTypeA/B 200:pre_postpre",
      "overallBlockTypeA/B 25:pre_postpre"), "df"]
  edf_conservative <- median(block_level_df)
  t_crit <- qt(0.975, edf_conservative)

  emmeans(fit.lmmPP, list(pairwise ~ overallBlockType), adjust = "tukey")

  emm_s.t <- emmeans(fit.lmmPP, pairwise ~ overallBlockType| mapStimLevel)
  emm_s.t <- emmeans(fit.lmmPP, pairwise ~ mapStimLevel | overallBlockType)
  summary(glht(fit.lmmPP,linfct=mcp(overallBlockType="Tukey")))

  emm_pairwise <- emmeans(fit.lmmPP,~overallBlockType*pre_post,adjust="Tukey")
  contrast(emm_pairwise,interaction="pairwise")
  # Effect sizes use total SD (RE variance + residual variance) so d values
  # are comparable to published Cohen's d benchmarks (small=0.2, medium=0.5, large=0.8).
  # Conservative df (block groups - fixed params) since df.residual returns
  # naive trial-level df (~6400).
  # Extract total variance from all RE components. With correlated random
  # slopes, VarCorr returns matrices, so sum diagonal (variances) not off-diagonal.
  vc <- VarCorr(fit.lmmPP)
  re_var <- sum(sapply(vc, function(v) sum(diag(v))))
  sigma_total <- sqrt(re_var + sigma(fit.lmmPP)^2)
  emm_effsize <- eff_size(emm_pairwise,sigma=sigma_total,edf=edf_conservative)

  # --- Plot: Interaction (overallBlockType x pre_post) with CIs ---
  # Use emmeans-native CIs (Satterthwaite df per cell, appropriate for
  # the cell-mean crossed RE model where each cell has its own df)
  emm_interaction <- emmeans(fit.lmmPP, ~ overallBlockType * pre_post)
  emm_int_df <- as.data.frame(emm_interaction)
  emm_int_df$pre_post <- factor(emm_int_df$pre_post, levels=c('pre','post'))
  pd <- position_dodge(width=0.15)
  marginal_means_plot <- ggplot(emm_int_df, aes(x=pre_post, y=emmean,
      color=overallBlockType, group=overallBlockType)) +
    geom_line(linewidth=1.2, position=pd) +
    geom_point(size=4, position=pd) +
    geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=0.06, linewidth=0.8, position=pd) +
    labs(x = "Pre versus Post Conditioning",
         y = expression(paste("Estimated Marginal Mean log(",mu,"V)")),
         color = "Experimental Condition",
         title = "Estimated Marginal Means by Conditioning and Status") +
    scale_color_brewer(palette="Dark2") +
    theme_bw(base_size=16) +
    theme(plot.title = element_text(size=16, face="bold"),
          axis.title = element_text(size=14),
          axis.text = element_text(size=12),
          legend.title = element_text(size=12),
          legend.text = element_text(size=11))
  print(marginal_means_plot)

  # --- Plot: EMM for each condition with CIs ---
  emm_condition <- emmeans(fit.lmmPP, ~ overallBlockType)
  emm_condition_df <- as.data.frame(emm_condition)
  p_emm_cond <- ggplot(emm_condition_df, aes(x=overallBlockType, y=emmean, color=overallBlockType)) +
    geom_point(size=3) +
    geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=0.2, linewidth=0.8) +
    labs(x = "Experimental Condition",
         y = expression(paste("Estimated Marginal Mean log(",mu,"V)")),
         color = "Condition",
         title = "Estimated Marginal Means by Conditioning Protocol") +
    scale_color_brewer(palette="Dark2") +
    theme_bw()
  print(p_emm_cond)

  # --- Plot: Pre-post contrast within each condition (95% CI) ---
  emm_prepost <- emmeans(fit.lmmPP, ~ pre_post | overallBlockType)
  prepost_contrasts <- contrast(emm_prepost, method="pairwise")
  prepost_df <- as.data.frame(confint(prepost_contrasts))
  p_prepost <- ggplot(prepost_df, aes(x=overallBlockType, y=estimate, color=overallBlockType)) +
    geom_point(size=3) +
    geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=0.2, linewidth=0.8) +
    geom_hline(yintercept=0, linetype="dashed", color="grey50") +
    labs(x = "Experimental Condition",
         y = expression(paste("Pre - Post Contrast log(",mu,"V)")),
         color = "Condition",
         title = "Pre vs Post Change by Conditioning Protocol (95% CI)") +
    scale_color_brewer(palette="Dark2") +
    theme_bw()
  print(p_prepost)

  # --- Plot: Effect sizes for pre-post contrasts ---
  prepost_effsize <- eff_size(emm_prepost, sigma=sigma_total, edf=edf_conservative)
  prepost_effsize_df <- as.data.frame(confint(prepost_effsize))
  p_effsize <- ggplot(prepost_effsize_df, aes(x=overallBlockType, y=effect.size, color=overallBlockType)) +
    geom_point(size=4) +
    geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=0.2, linewidth=1.0) +
    geom_hline(yintercept=0, linetype="dashed", color="grey50") +
    labs(x = "Experimental Condition",
         y = "Standardized Effect Size (d)",
         color = "Condition",
         title = "Effect Sizes: Pre vs Post by Conditioning Protocol (95% CI)",
         subtitle = "Computed on log-transformed EP magnitude") +
    scale_color_brewer(palette="Dark2") +
    theme_bw(base_size=16) +
    theme(plot.title = element_text(size=16, face="bold"),
          plot.subtitle = element_text(size=12),
          axis.title = element_text(size=14),
          axis.text = element_text(size=12),
          legend.title = element_text(size=12),
          legend.text = element_text(size=11))
  print(p_effsize)

  # --- Plot: Compare fixed effects across all three model approaches ---
  extract_coefs <- function(fit, label) {
    df <- as.data.frame(summary(fit)$coefficients)
    df$term <- rownames(df)
    df$model <- label
    df
  }
  coef_both <- rbind(
    extract_coefs(fit.lmmPP, "1: Channel slope (primary)"),
    extract_coefs(fit.lmmPP.blockslope, "2: Block slope"),
    extract_coefs(fit.lmmPP.intonly, "3: Intercept-only"),
    extract_coefs(fit.lmmPP.trial.flat, "4: Trial-level flat RE")
  )
  coef_both <- coef_both %>% filter(term != "(Intercept)")
  p_compare <- ggplot(coef_both, aes(x=term, y=Estimate, color=model)) +
    geom_point(position=position_dodge(width=0.6), size=2.5) +
    geom_errorbar(aes(ymin=Estimate - t_crit*`Std. Error`, ymax=Estimate + t_crit*`Std. Error`),
                  position=position_dodge(width=0.6), width=0.3) +
    geom_hline(yintercept=0, linetype="dashed", color="grey50") +
    labs(x = "", y = "Estimate (95% CI)", color = "Model",
         title = "Fixed Effect Estimates: All Three Approaches") +
    theme_bw() +
    coord_flip()
  print(p_compare)

  emmFigSize = 7  # square dimensions for EMM plots
  if (savePlot && !avgMeas) {
    ggsave(plot=marginal_means_plot, paste0("emmip_AUC.png"), units="in", width=emmFigSize, height=emmFigSize, dpi=300)
    ggsave(plot=marginal_means_plot, paste0("emmip_AUC.eps"), units="in", width=emmFigSize, height=emmFigSize, device=cairo_ps, fallback_resolution=300)
    ggsave(plot=p_emm_cond, paste0("emm_condition_AUC.png"), units="in", width=emmFigSize, height=emmFigSize, dpi=300)
    ggsave(plot=p_emm_cond, paste0("emm_condition_AUC.eps"), units="in", width=emmFigSize, height=emmFigSize, device=cairo_ps, fallback_resolution=300)
    ggsave(plot=p_prepost, paste0("prepost_contrast_AUC.png"), units="in", width=emmFigSize, height=emmFigSize, dpi=300)
    ggsave(plot=p_prepost, paste0("prepost_contrast_AUC.eps"), units="in", width=emmFigSize, height=emmFigSize, device=cairo_ps, fallback_resolution=300)
    ggsave(plot=p_effsize, paste0("prepost_effsize_AUC.png"), units="in", width=emmFigSize, height=emmFigSize, dpi=300)
    ggsave(plot=p_effsize, paste0("prepost_effsize_AUC.eps"), units="in", width=emmFigSize, height=emmFigSize, device=cairo_ps, fallback_resolution=300)
    ggsave(plot=p_compare, paste0("model_comparison_AUC.png"), units="in", width=figWidth, height=figHeight, dpi=300)
    ggsave(plot=p_compare, paste0("model_comparison_AUC.eps"), units="in", width=figWidth, height=figHeight, device=cairo_ps, fallback_resolution=300)

  } else if (savePlot && avgMeas){
    ggsave(plot=marginal_means_plot, paste0("emmip_avg_AUC.png"), units="in", width=emmFigSize, height=emmFigSize, dpi=300)
    ggsave(plot=marginal_means_plot, paste0("emmip_avg_AUC.eps"), units="in", width=emmFigSize, height=emmFigSize, device=cairo_ps, fallback_resolution=300)

  }
  
  }
  else {
  fit.lmmPP = lmerTest::lmer(PPvec ~ mapStimLevel + disease  + chanInCond + blockType + baLabel + (1|subjectNum/chanVec),data=dataListAgg)
  emmeans(fit.lmmPP, list(pairwise ~ blockType), adjust = "tukey")
  
  emm_s.t <- emmeans(fit.lmmPP, pairwise ~ blockType | mapStimLevel)
  emm_s.t <- emmeans(fit.lmmPP, pairwise ~ mapStimLevel | blockType)
  summary(glht(fit.lmmPP,linfct=mcp(blockType="Tukey")))
  
  }
  ## this one is for only the highest stimulation level 
  # fit.lmmPP = lmerTest::lmer(PPvec ~ disease  + chanInCond + blockType + (1|subjectNum),data=dataList)
  
  # fit.lmmPP = lm(PPvec ~ mapStimLevel + disease  + chanInCond + blockType ,data=dataList)
  
  #fit.lmmPP = lmerTest::lmer(PPvec ~ mapStimLevel + blockType + disease + (1|sidVec:chanInCond),data=dataList)
  #fit.lmmPP = lmerTest::lmer(PPvec ~ mapStimLevel + blockType + disease + (1|subjectNum),data=dataList)
  
#emm_options(pbkrtest.limit = 200000) 
  
  summary(fit.lmmPP)
  plot(fit.lmmPP)
  qqnorm(resid(fit.lmmPP))
  qqline(resid(fit.lmmPP))
  summary(glht(fit.lmmPP,linfct=mcp(mapStimLevel="Tukey")))

  #emm_options(pbkrtest.limit = 10000)

  emmeans(fit.lmmPP, list(pairwise ~ mapStimLevel), adjust = "tukey")
  #emmeans(fit.lmmPP, list(pairwise ~ baLabel), adjust = "tukey") # baLabel not in reduced model


  cat("\n=== TYPE III ANOVA (primary model) ===\n")
  print(anova(fit.lmmPP, type=3))
  cat("\n=== PRE-POST CONTRASTS ===\n")
  emm_pp2 <- emmeans(fit.lmmPP, ~ pre_post | overallBlockType)
  print(confint(contrast(emm_pp2, method="pairwise")))
  cat("\n=== TUKEY PAIRWISE (conditioning protocol) ===\n")
  emm_cond2 <- emmeans(fit.lmmPP, ~ overallBlockType)
  print(pairs(emm_cond2, adjust="tukey"))
  # Residual diagnostics
  cat("\n=== RESIDUAL DIAGNOSTICS (primary model) ===\n")
  r <- resid(fit.lmmPP)
  cat("N residuals:", length(r), "\n")
  sw <- shapiro.test(r)
  cat("Shapiro-Wilk W =", round(sw$statistic, 4), ", p =", format(sw$p.value, digits=3), "\n")
  n <- length(r)
  skew <- (sum((r - mean(r))^3) / n) / (sum((r - mean(r))^2) / n)^1.5
  kurt <- (sum((r - mean(r))^4) / n) / (sum((r - mean(r))^2) / n)^2 - 3
  cat("Skewness:", round(skew, 3), "\n")
  cat("Excess kurtosis:", round(kurt, 3), "\n")

  cat("\n=== EFFECT SIZES ===\n")
  cat("Total SD:", sigma_total, "\n")
  es2 <- eff_size(emm_pp2, sigma=sigma_total, edf=edf_conservative)
  print(confint(es2))
  tab_model(fit.lmmPP)

  # --- Exploratory: baLabel as fixed effect ---
  cat("\n=== EXPLORATORY: baLabel distribution ===\n")
  print(table(dataListAgg$baLabel))
  ba_per_chan <- dataListAgg %>% group_by(subjectNum, chanVec) %>%
    dplyr::summarise(ba=first(as.character(baLabel)), .groups="drop")
  print(ba_per_chan)

  cat("\n=== EXPLORATORY: baLabel as fixed effect ===\n")
  fit.ba <- tryCatch({
    lmerTest::lmer(PPvec ~ mapStimLevel + chanInCond + baLabel +
      overallBlockType*pre_post + (1|subjectNum/chanVec) + (1|subjectNum:blockVec),
      data=dataListAgg)
  }, error=function(e) { cat("ERROR:", e$message, "\n"); NULL })
  if (!is.null(fit.ba)) {
    cat("Singular:", isSingular(fit.ba), "\n")
    print(summary(fit.ba))
    cat("\n=== ANOVA with baLabel ===\n")
    print(anova(fit.ba, type=3))
  }

  # Random slope model is now the primary model (see fit.lmmPP above).
  # LRT vs intercept-only printed in the comparison section.

  # Channel slope is primary. Block slope and both-slopes tested during
  # model selection: block slope AIC=210, channel slope AIC=158, both
  # slopes AIC=161 (block slope collapses to SD~0.04 when channel slope
  # is present). Dose-response steepness is a property of electrode
  # location, not recording session.

  # Save side-by-side comparison of all model approaches
  tab_comparison <- tab_model(fit.lmmPP, fit.lmmPP.blockslope, fit.lmmPP.intonly, fit.lmmPP.trial.flat,
    show.re.var=TRUE, show.icc=TRUE, show.obs=TRUE,
    dv.labels=c("Channel slope (primary)","Block slope","Intercept-only","Trial-level flat RE"))
  writeLines(as.character(tab_comparison$page.complete), paste0(outputDir, "/tab_model_comparison.html"))

  # Also save tab_model as Word-compatible HTML (can open directly in Word)
  tab_primary <- tab_model(fit.lmmPP, show.re.var=TRUE, show.icc=TRUE, show.obs=TRUE,
    dv.labels="Primary Model (crossed RE + channel dose slope)",
    title="Linear Mixed-Effects Model Results")
  writeLines(as.character(tab_primary$page.complete), paste0(outputDir, "/tab_model_primary.html"))
  cat("Saved tab_model HTML (Word-compatible) to:", paste0(outputDir, "/tab_model_primary.html"), "\n")

  # --- Generate manuscript-ready .docx tables ---
  if (requireNamespace("officer", quietly=TRUE) && requireNamespace("flextable", quietly=TRUE)) {
    library(officer)
    library(flextable)

    fmt_p <- function(p) ifelse(p < 0.001, "< 0.001", sprintf("%.3f", p))

    doc <- read_docx()

    # Table: Random Effects
    vc <- VarCorr(fit.lmmPP)
    ngrps <- summary(fit.lmmPP)$ngrps
    re_rows <- list()
    for (nm in names(vc)) {
      v <- vc[[nm]]
      ng <- as.character(ngrps[nm])
      if (ncol(v) == 1) {
        re_rows[[length(re_rows)+1]] <- data.frame(Component=nm, Term="(Intercept)",
          Variance=round(v[1,1],4), SD=round(sqrt(v[1,1]),4), Corr="", Groups=ng, stringsAsFactors=FALSE)
      } else {
        corr_val <- sprintf("%.2f", attr(v,"correlation")[2,1])
        re_rows[[length(re_rows)+1]] <- data.frame(
          Component=c(nm,""), Term=rownames(v),
          Variance=round(diag(v),4), SD=round(sqrt(diag(v)),4),
          Corr=c("", corr_val), Groups=c(ng,""), stringsAsFactors=FALSE)
      }
    }
    re_df <- do.call(rbind, re_rows)
    re_df <- rbind(re_df, data.frame(Component="Residual", Term="", Variance=round(sigma(fit.lmmPP)^2,4),
                                      SD=round(sigma(fit.lmmPP),4), Corr="", Groups=""))
    doc <- body_add_par(doc, "Table: Random Effects", style="heading 2")
    ft <- flextable(re_df) |> autofit() |>
      set_caption("Random effects from primary model (cell-median, crossed RE + random linear dose slope)")
    doc <- body_add_flextable(doc, ft)
    doc <- body_add_par(doc, sprintf("Model is not singular. AIC = %.1f", AIC(fit.lmmPP)))
    doc <- body_add_par(doc, "")

    # Table: Type III ANOVA
    aov_tbl <- as.data.frame(anova(fit.lmmPP, type=3))
    aov_tbl$Effect <- rownames(aov_tbl)
    aov_tbl <- aov_tbl[, c("Effect","Sum Sq","Mean Sq","NumDF","DenDF","F value","Pr(>F)")]
    aov_tbl$`Sum Sq` <- round(aov_tbl$`Sum Sq`, 3)
    aov_tbl$`Mean Sq` <- round(aov_tbl$`Mean Sq`, 3)
    aov_tbl$DenDF <- round(aov_tbl$DenDF, 1)
    aov_tbl$`F value` <- round(aov_tbl$`F value`, 2)
    aov_tbl$p <- sapply(aov_tbl$`Pr(>F)`, fmt_p)
    aov_tbl$`Pr(>F)` <- NULL
    doc <- body_add_par(doc, "Table: Type III ANOVA (Satterthwaite df)", style="heading 2")
    ft <- flextable(aov_tbl) |> autofit() |>
      set_caption("Type III ANOVA with Satterthwaite denominator df")
    doc <- body_add_flextable(doc, ft)
    doc <- body_add_par(doc, "")

    # Table: Fixed Effects
    fe_tbl <- as.data.frame(summary(fit.lmmPP)$coefficients)
    fe_tbl$Predictor <- rownames(fe_tbl)
    fe_tbl <- fe_tbl[, c("Predictor","Estimate","Std. Error","df","t value","Pr(>|t|)")]
    fe_tbl$Estimate <- round(fe_tbl$Estimate, 3)
    fe_tbl$`Std. Error` <- round(fe_tbl$`Std. Error`, 3)
    fe_tbl$df <- round(fe_tbl$df, 1)
    fe_tbl$`t value` <- round(fe_tbl$`t value`, 2)
    fe_tbl$p <- sapply(fe_tbl$`Pr(>|t|)`, fmt_p)
    fe_tbl$`Pr(>|t|)` <- NULL
    doc <- body_add_par(doc, "Table: Fixed Effects", style="heading 2")
    ft <- flextable(fe_tbl) |> autofit() |>
      set_caption("Fixed effects. All estimates on log(uV) scale.")
    doc <- body_add_flextable(doc, ft)
    doc <- body_add_par(doc, "")

    # Table: Pre-Post Contrasts
    emm_pp_docx <- emmeans(fit.lmmPP, ~ pre_post | overallBlockType)
    contrast_tbl <- as.data.frame(confint(contrast(emm_pp_docx, method="pairwise")))
    es_tbl <- as.data.frame(confint(eff_size(emm_pp_docx, sigma=sigma_total, edf=edf_conservative)))
    contrast_tbl$`Cohen's d` <- round(es_tbl$effect.size, 2)
    contrast_tbl$`d lower` <- round(es_tbl$lower.CL, 3)
    contrast_tbl$`d upper` <- round(es_tbl$upper.CL, 3)
    contrast_tbl$`d df` <- round(es_tbl$df, 1)
    contrast_tbl$estimate <- round(contrast_tbl$estimate, 3)
    contrast_tbl$SE <- round(contrast_tbl$SE, 3)
    contrast_tbl$df <- round(contrast_tbl$df, 1)
    contrast_tbl$lower.CL <- round(contrast_tbl$lower.CL, 3)
    contrast_tbl$upper.CL <- round(contrast_tbl$upper.CL, 3)
    contrast_out <- contrast_tbl[, c("overallBlockType","estimate","SE","lower.CL","upper.CL","df","Cohen's d","d lower","d upper","d df")]
    names(contrast_out)[1] <- "Condition"
    names(contrast_out)[4:5] <- c("CI lower", "CI upper")
    doc <- body_add_par(doc, "Table: Pre-Post Contrasts", style="heading 2")
    ft <- flextable(contrast_out) |> autofit() |>
      set_caption(sprintf("EMM contrasts (post - pre). Cohen's d uses total SD = %.3f.", sigma_total))
    doc <- body_add_flextable(doc, ft)
    doc <- body_add_par(doc, "")

    # Table: Tukey Pairwise
    emm_cond_docx <- emmeans(fit.lmmPP, ~ overallBlockType)
    tukey_tbl <- as.data.frame(pairs(emm_cond_docx, adjust="tukey"))
    tukey_tbl$estimate <- round(tukey_tbl$estimate, 3)
    tukey_tbl$SE <- round(tukey_tbl$SE, 3)
    tukey_tbl$df <- round(tukey_tbl$df, 1)
    tukey_tbl$t.ratio <- round(tukey_tbl$t.ratio, 2)
    tukey_tbl$p.value <- sapply(tukey_tbl$p.value, fmt_p)
    names(tukey_tbl)[1] <- "Comparison"
    doc <- body_add_par(doc, "Table: Pairwise Comparisons (Tukey)", style="heading 2")
    ft <- flextable(tukey_tbl) |> autofit() |>
      set_caption("Pairwise comparisons of conditioning protocols (Tukey-adjusted)")
    doc <- body_add_flextable(doc, ft)

    docx_path <- paste0(outputDir, "/statistical_tables.docx")
    print(doc, target=docx_path)
    cat("Saved manuscript tables to:", docx_path, "\n")
  }

  if (savePlot && !avgMeas) {
    
    figHeight = 4
    figWidth = 8
    png("pairedPulse_resid_PP_trim_allSubjs_AUC.png",width=figWidth,height=figHeight,units="in",res=600)
    plot(fit.lmmPP)
    dev.off()

    setEPS()
    postscript("pairedPulse_resid_PP_trim_allSubjs_AUC.eps",width=figWidth,height=figHeight)
    plot(fit.lmmPP)
    dev.off()
    
    figHeight = 4
    figWidth = 8
    png("pairedPulse_qq_PP__trim_allSubjs_AUC.png",width=figWidth,height=figHeight,units="in",res=600)
    qqnorm(resid(fit.lmmPP))
    qqline(resid(fit.lmmPP))  #summary(fit.lmm2)
    dev.off()
    
    setEPS()
    postscript("pairedPulse_qq_PP_trim_allSubjs_AUC.eps",width=figWidth,height=figHeight)
    qqnorm(resid(fit.lmmPP))
    qqline(resid(fit.lmmPP))  #summary(fit.lmm2)
    dev.off()
    
  } else if (savePlot && avgMeas){
    
    figHeight = 4
    figWidth = 8
    png("pairedPulse_resid_PP_allSubjs_avg_AUC.png",width=figWidth,height=figHeight,units="in",res=600)
    plot(fit.lmmPP)
    dev.off()

    setEPS()
    postscript("pairedPulse_resid_PP_allSubjs_avg_AUC.eps",width=figWidth,height=figHeight)
    plot(fit.lmmPP)
    dev.off()
    
    figHeight = 4
    figWidth = 8
    png("pairedPulse_qq_PP_allSubjs_avg_AUC.png",width=figWidth,height=figHeight,units="in",res=600)
    qqnorm(resid(fit.lmmPP))
    qqline(resid(fit.lmmPP))  #summary(fit.lmm2)
    dev.off()
    
    setEPS()
    postscript("pairedPulse_qq_PP_allSubjs_avg_AUC.eps",width=figWidth,height=figHeight)
    qqnorm(resid(fit.lmmPP))
    qqline(resid(fit.lmmPP))  #summary(fit.lmm2)
    dev.off()
    
  }
  
  # Reset emmeans df method for small-data models (Satterthwaite is fast on 120 rows)
  emm_options(lmer.df = "satterthwaite")

  #fit.lmmdiff = lmerTest::lmer(absDiff ~ stimLevelVec + blockType + chanInCond + (1|sidVec),data=dataList)
  fit.lmmdiff = lmerTest::lmer(absDiff ~ mapStimLevel + chanInCond + blockType + (1|subjectNum/chanVec),data=dataListAgg)


  fit.lm = lm(absDiff ~ mapStimLevel + blockType + chanInCond ,data=dataListAgg)
  
  summary(fit.lmmdiff)
  plot(fit.lmmdiff)
  
  qqnorm(resid(fit.lmmdiff))
  qqline(resid(fit.lmmdiff))  #summary(fit.lmm2)
  summary(glht(fit.lmmdiff,linfct=mcp(blockType="Tukey")))
  summary(glht(fit.lmmdiff,linfct=mcp(mapStimLevel="Tukey")))
  
  emmeans(fit.lmmdiff, list(pairwise ~ blockType), adjust = "tukey")
  emmeans(fit.lmmdiff, list(pairwise ~ mapStimLevel), adjust = "tukey")
  
  emm_s.t <- emmeans(fit.lmmdiff, pairwise ~ blockType | mapStimLevel)
  emm_s.t <- emmeans(fit.lmmdiff, pairwise ~ mapStimLevel | blockType)
  
  anova(fit.lmmdiff)
  tab_model(fit.lmmdiff)
  
  if (savePlot && !avgMeas) {
    
    figHeight = 4
    figWidth = 8
    png("pairedPulse_resid_diff_allSubjs_AUC.png",width=figWidth,height=figHeight,units="in",res=600)
    plot(fit.lmmdiff)
    dev.off()
    
    setEPS()
    postscript("pairedPulse_resid_diff_allSubjs_AUC.eps",width=figWidth,height=figHeight)
    plot(fit.lmmdiff)
    dev.off()
    
    figHeight = 4
    figWidth = 8
    png("pairedPulse_qq_diff_allSubjs_AUC.png",width=figWidth,height=figHeight,units="in",res=600)
    qqnorm(resid(fit.lmmdiff))
    qqline(resid(fit.lmmdiff))  #summary(fit.lmm2)
    dev.off()
    
    setEPS()
    postscript("pairedPulse_qq_diff_allSubjs_AUC.eps",width=figWidth,height=figHeight)
    qqnorm(resid(fit.lmmdiff))
    qqline(resid(fit.lmmdiff))  #summary(fit.lmm2)
    dev.off()
    
  } else if (savePlot && avgMeas){
    
    figHeight = 4
    figWidth = 8
    png("pairedPulse_qq_resid_diff_allSubjs_avg_AUC.png",width=figWidth,height=figHeight,units="in",res=600)
    plot(fit.lmmdiff)
    dev.off()
    
    setEPS()
    postscript("pairedPulse_resid_diff_allSubjs_avg_AUC.eps",width=figWidth,height=figHeight)
    plot(fit.lmmdiff)
    dev.off()
    
    figHeight = 4
    figWidth = 8
    png("pairedPulse_qq_diff_allSubjs_avg_AUC.png",width=figWidth,height=figHeight,units="in",res=600)
    qqnorm(resid(fit.lmmdiff))
    qqline(resid(fit.lmmdiff))  #summary(fit.lmm2)
    dev.off()
    
    setEPS()
    postscript("pairedPulse_qq_diff_allSubjs_avg_AUC.eps",width=figWidth,height=figHeight)
    qqnorm(resid(fit.lmmdiff))
    qqline(resid(fit.lmmdiff))  #summary(fit.lmm2)
    dev.off()
    
  }
  
  # Per-cell trial-level Cohen's d analysis removed. With ~5-10 trials per
  # cell on log-transformed data, per-cell d values are extremely noisy
  # (d values up to +-10 from near-zero denominators). Model-based effect
  # sizes via emmeans::eff_size() with total SD denominator (prepost_effsize_AUC)
  # are the principled approach and are reported above.

}

# Save workspace for quick loading without rerunning
save.image(file=here("DBS_EP_PairedPulse","R_output","main_analysis_workspace.RData"))
cat("Saved workspace to R_output/main_analysis_workspace.RData\n")

setwd(oldwd)
