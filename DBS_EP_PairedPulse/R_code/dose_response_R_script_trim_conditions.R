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
library('effectsize')
library('caret')
library('DescTools')

rootDir = here()
dataDir = here("DBS_EP_PairedPulse","R_data")
codeDir = here("DBS_EP_PairedPulse","R_code")

sidVec <- c('46c2a','c963f','2e114','fe7df','e6f3c','9f852',
            '8e907','08b13','e9c9b','41a73','68574',
            '01fee','a23ed')

diseaseVec <- c('PD','PD','MD','PD','PD','PD','PD','MD',
                'PD','PD','PD','PD','MD')

sidVec <- c('c963f','2e114','fe7df','e6f3c',
            '8e907','08b13','e9c9b','41a73','68574',
            '01fee')

diseaseVec <- c('PD','MD','PD','PD','PD','MD',
                'PD','PD','PD','PD')



sidVec <- c('46c2a','c963f','2e114','fe7df','e6f3c','9f852',
            '8e907','08b13','e9c9b','41a73','68574',
            '01fee')


diseaseVec <- c('MD','PD','MD','PD','PD','PD','PD','MD',
                'PD','PD','PD','PD')

#sidVec <- c('9f852')
#sidVec <- c('41a73')
#diseaseVec <- c('PD')
repeatedMeasures = TRUE # if true, does repeated measures analysis, if false, does more of ANCOVA style analysis
log_data = TRUE
box_data = FALSE
trim_data = TRUE
min_stim_level = 3
savePlot = 0
avgMeasVec = c(0)
figWidth = 8 
figHeight = 6 

for (avgMeas in avgMeasVec) {
  
  dataList = list()
  blockList = list()
  conditionList = list()
  index = 1
  
  subjectNumIndex = 1
  
  for (sid in sidVec){
    source(here("DBS_EP_PairedPulse","R_config_files",paste0("subj_",sid,".R")))
    
    brodmann_areas <- read.csv(here("DBS_EP_PairedPulse","R_config_files",paste0(sid,'_MNIcoords_labelled.csv')),header=TRUE,sep = ",",stringsAsFactors=F)
    if (avgMeas) {
      dataPP <- read.table(here("DBS_EP_PairedPulse","R_data",paste0(sid,'_PairedPulseData_new_pk_pk_avg_5.csv')),header=TRUE,sep = ",",stringsAsFactors=F, colClasses=c("stimLevelVec"="numeric","sidVec"="factor"))
    } else{
      dataPP <- read.table(here("DBS_EP_PairedPulse","R_data",paste0(sid,'_PairedPulseData_new_rms_pk_pk.csv')),header=TRUE,sep = ",",stringsAsFactors=F, colClasses=c("stimLevelVec"="numeric","sidVec"="factor"))
    }
    
    if (sid == "2e114"){
      goodVec = c(1000,2500,3500)
      dataPP <- dataPP %>% filter(stimLevelVec %in% goodVec)
    }
    
    # multiply by 1e6
    # here is where we use the column from the average one
    dataPP$PPvec = dataPP$rmsVec*1e6
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
      
      dataInt$baLabel = as.factor(dataInt$baLabel)
      dataInt$aalLabel = as.factor(dataInt$aalLabel)
 
     # dataInt <- subset(dataInt, PPvec<1000)
      #dataInt <- subset(dataInt, PPvec>30)
      
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
      

      for (comparison in whichCompareVec){
        dataIntCompare <- subset(dataInt,(blockVec %in% comparison))

        #if(!(is.data.frame(dataIntCompare) & (nrow(dataIntCompare)==0))){

        # make sure each part of the comparison has at least 10 trials, otherwise exclude it 
        if (avgMeas){
          dataIntCompare <- dataIntCompare %>% group_by(mapStimLevel,blockVec) %>% mutate(enough = n()>=2)
          
        }
        else {
          dataIntCompare <- dataIntCompare %>% group_by(mapStimLevel,blockVec) %>% mutate(enough = n()>=10)
          
        }
        dataIntCompare <- dataIntCompare %>% group_by(mapStimLevel) %>% filter(all(enough) & length(unique(blockVec))>1)
        
        if(!(empty(dataIntCompare))){   
          
        dataIntCompare <- dataIntCompare %>% 
          group_by(mapStimLevel) %>% 
          mutate(
            effectSize = cohens_d(PPvec[blockVec==comparison[2]],PPvec[blockVec==comparison[1]])[1,1],
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
      
      p <- ggplot(dataInt, aes(x=stimLevelVec, y=PPvec,color=stimLevelVec)) +
        geom_point(position=position_jitterdodge(dodge.width=0.250)) +  geom_smooth(method=lm) + facet_wrap(~blockVec, labeller = as_labeller(blockNames))+
        labs(x = expression(paste("Stimulation Current (mA)")),y=expression(paste("Peak to Peak Magnitude (",mu,"V)")),title = paste0("Subject ", subjectNumInterest," DBS Paired Pulse EP Measurements")) +
        guides(colour=guide_colorbar("Stimulation Level"))
      print(p)
      
      if (savePlot && !avgMeas) {
        ggsave(paste0("subj_", subjectNum, "_ID_", sid,"_chan_",chanInt,"_scatter_lm.png"), units="in", width=figWidth, height=figHeight, dpi=600)
        ggsave(paste0("subj_", subjectNum, "_ID_", sid,"_chan_",chanInt,"_scatter_lm.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
        
      }
      else if (savePlot && avgMeas){
        ggsave(paste0("subj_", subjectNum, "_ID_", sid,"_chan_",chanInt,"_scatter_lm_avg.png"), units="in", width=figWidth, height=figHeight, dpi=600)
        ggsave(paste0("subj_", subjectNum, "_ID_", sid,"_chan_",chanInt,"_scatter_lm_avg.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
        
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
    subjectNumIndex = subjectNumIndex + 1
  }
  
  # now do comparison at highest stim level relative to baseline
  dataList <- do.call(rbind,dataList)
  #blockList <- do.call(rbind,blockList)
  
  if (repeatedMeasures){
    dataList <- dataList %>% filter((blockType %in% c('baseline','A/B 25','A/B 200','A/A 200')) & (overallBlockType %in% c('A/B 25','A/B 200','A/A 200')))
  }
  else {
    dataList <- dataList %>% filter(blockType %in% c('baseline','A/B 25','A/B 200','A/A 200'))
  }
  #plot
  grouped <- group_by(dataList, sidVec, chanVec, blockType,mapStimLevel,disease)
  dataListSummarize <- summarise(grouped,meanPerc = mean(percentDiff),sdPerc = sd(percentDiff),
                                 meanAbs=mean(absDiff), sdDiff=sd(percentDiff),meanPP = mean(PPvec),sdPP = sd(PPvec))
  
  p3 <- ggplot(dataList, aes(x=as.numeric(mapStimLevel), y=percentDiff,color=blockType)) +
    geom_point(position=position_jitterdodge(dodge.width=0.75)) +geom_smooth(method=lm) +
    labs(x = expression(paste("Stimulation Level")),y=expression(paste("Percent Difference in EP Magnitude from Baseline")), color="Experimental Condition",title = paste0("Changes in EP Magnitude"))
  print(p3)
  
  if (savePlot && !avgMeas) {
    ggsave(paste0("across_subj_percent.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    ggsave(paste0("across_subj_percent.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
    
  } else if (savePlot && avgMeas){
    ggsave(paste0("across_subj_percent_avg.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    ggsave(paste0("across_subj_percent_avg.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
    
  }
  
  p4 <- ggplot(dataList, aes(x=as.numeric(mapStimLevel), y=absDiff,color=blockType)) +
    geom_point(position=position_jitterdodge(dodge.width=0.75)) +geom_smooth(method=lm) +
    labs(x = expression(paste("Stimulation Level")),y=expression(paste("Absolute Difference in EP Magnitude from Baseline (",mu,"V)")),color="Experimental Condition",title = paste0("Changes in EP Magnitude"))
  print(p4)
  
  if (savePlot && !avgMeas) {
    ggsave(paste0("across_subj_abs.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    ggsave(paste0("across_subj_abs.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
    
  } else if (savePlot && avgMeas){
    ggsave(paste0("across_subj_abs_avg.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    ggsave(paste0("across_subj_abs_avg.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
    
  }
  
  p5 <- ggplot(dataListSummarize, aes(x=as.numeric(mapStimLevel), y=meanPerc,color=blockType)) +
    geom_point(position=position_jitterdodge(dodge.width=0.5)) +geom_smooth(method=lm) +
    labs(x = expression(paste("Stimulation Level")),y=expression(paste("Mean Percent Difference in EP Magnitude")), color="Experimental Condition",title = paste0("Changes in EP Magnitude"))
  print(p5)
  
  if (savePlot && !avgMeas) {
    ggsave(paste0("across_subj_mean_perc.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    ggsave(paste0("across_subj_mean_perc.eps"), units="in", width=figWidth, height=figHeight, device=cairo_ps, fallback_resolution=600)
    
  } else if (savePlot && avgMeas){
    ggsave(paste0("across_subj_mean_perc_avg.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    ggsave(paste0("across_subj_mean_perc_avg.eps"), units="in", width=figWidth, height=figHeight, device=cairo_ps, fallback_resolution=600)
    
  }
  
  p6 <- ggplot(dataListSummarize, aes(x=mapStimLevel, y=meanAbs,color=blockType)) +
    geom_point(position=position_jitterdodge(dodge.width=0.5)) +geom_smooth(method=lm) +
    labs(x = expression(paste("Stimulation Level")),y=expression(paste("Mean Absolute Difference in EP Magnitude from Baseline (",mu,"V)")),
         color="Experimental Condition",title = paste0("Changes in EP Magnitude"))
  print(p6)
  
  if (savePlot && !avgMeas) {
    ggsave(paste0("across_subj_mean_abs.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    ggsave(paste0("across_subj_mean_abs.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
    
  } else if (savePlot && avgMeas){
    ggsave(paste0("across_subj_mean_abs_avg.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    ggsave(paste0("across_subj_mean_abs_avg.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
    
  }
  
  p7 <- ggplot(dataList, aes(x=mapStimLevel, y=PPvec,color=blockType)) +
    geom_point(position=position_jitterdodge(dodge.width=0.75)) +geom_smooth(method=lm) +
    labs(x = expression(paste("Stimulation Level")),y=expression(paste("Peak to Peak Voltage (",mu,"V)")),
         color="Experimental Condition",title = paste0("EP Magnitude"))
  print(p7)
  if (savePlot && !avgMeas) {
    ggsave(paste0("across_subj_PP.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    ggsave(paste0("across_subj_PP.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
    
  } else if (savePlot && avgMeas){
    ggsave(paste0("across_subj_PP_avg.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    ggsave(paste0("across_subj_PP_avg.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
    
  }

  p8 <- ggplot(dataListSummarize, aes(x=mapStimLevel, y=meanPP,color=blockType)) +
    geom_point(position=position_jitterdodge(dodge.width=0.5)) +geom_smooth(method=lm) +
    labs(x = expression(paste("Stimulation Level")),y=expression(paste("Mean peak to peak voltage (",mu,"V)")), color="Experimental condition",title = paste0("EP magnitude"))
  print(p8)
  
  if (savePlot && !avgMeas) {
    ggsave(paste0("across_subj_mean_PP.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    ggsave(paste0("across_subj_mean_PP.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
    
  } else if (savePlot && avgMeas){
    ggsave(paste0("across_subj_mean_PP_avg.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    ggsave(paste0("across_subj_mean_PP_avg.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
    
  }  
  
  # figWidth = 8
  # figHeight = 4
  # p11 <- ggplot(data = dataListSummarize, aes(x = blockType, y = meanAbs,color=blockType)) +
  #   geom_boxplot(notch=TRUE,outlier.shape=NA)  + geom_jitter(shape=16, position=position_jitter(0.2),aes(alpha = mapStimLevel)) +
  #   labs(x = expression(paste("Experimental Condition")),y=expression(paste("Absolute Difference from Baseline Peak-To-Peak (",mu,"V)")),color="Experimental Condition",alpha="Ordered Stim Level",title = paste0("EP Difference from Baseline by Conditioning Protocol"))
  # print(p11)
  # 
  # if (savePlot && !avgMeas) {
  #   ggsave(paste0("across_subj_mean_abs_box.png"), units="in", width=figWidth, height=figHeight, dpi=600)
  #   ggsave(paste0("across_subj_mean_abs_box.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
  # 
  # } else if (savePlot && avgMeas){
  #   ggsave(paste0("across_subj_mean_abs_box_avg.png"), units="in", width=figWidth, height=figHeight, dpi=600)
  #   ggsave(paste0("across_subj_mean_abs_box_avg.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
  # 
  # }
  
  # box plot
  
  figWidth = 8
  figHeight = 4
  goodVecBlock <- c("A/B 25","A/B 200","A/A 200")
  p11 <- ggplot(data =  dataListSummarize%>% filter(blockType %in% goodVecBlock), aes(x = blockType, y = meanAbs,color=blockType)) +
    geom_boxplot(notch=TRUE,outlier.shape=NA)  + geom_jitter(shape=16, position=position_jitter(0.2),aes(alpha = mapStimLevel)) +
    labs(x = expression(paste("Experimental Condition")),y=expression(paste("Absolute Difference from Baseline Peak-To-Peak (",mu,"V)")),color="Experimental Condition",alpha="Ordered Stim Level",title = paste0("EP Difference from Baseline by Conditioning Protocol")) +
    scale_color_brewer(palette="Dark2")
  print(p11)
  
  if (savePlot && !avgMeas) {
    ggsave(paste0("across_subj_mean_abs_box_all_subjs.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    ggsave(paste0("across_subj_mean_abs_box_all_subjs.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
    
  } else if (savePlot && avgMeas){
    ggsave(paste0("across_subj_mean_abs_box_avg.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    ggsave(paste0("across_subj_mean_abs_box_avg.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
    
  }
  dataBarInt = dataList%>% filter(blockType %in% goodVecBlock)
  groupedAll <- group_by(dataBarInt,blockType)
  dataListSummarizeBar <- summarise(groupedAll,meanPerc = mean(percentDiff),sdPerc = sd(percentDiff),
                                    meanAbs=mean(absDiff),sdAbs=sd(absDiff), sdDiff=sd(percentDiff),meanPP = mean(PPvec),sdPP = sd(PPvec),
                                    n=n(),seAbs=sdAbs/sqrt(n))
  
  # bar plot with error
  figWidth = 8
  figHeight = 4
  goodVecBlock <- c("A/B 25","A/B 200","A/A 200")
  p12 <- ggplot(data =  dataListSummarizeBar, aes(x = blockType, y = meanAbs,fill=blockType)) +
    geom_bar(stat="identity")  +
    labs(x = expression(paste("Experimental Condition")),y=expression(paste("Absolute Difference from Baseline Peak-To-Peak (",mu,"V)")),fill="Experimental Condition",title = paste0("EP Difference from Baseline by Conditioning Protocol")) +
    scale_fill_brewer(palette="Dark2") +
    geom_errorbar(aes(ymin=meanAbs-seAbs, ymax=meanAbs+seAbs), width=.2,
                  position=position_dodge(.9)) 
  print(p12)
  
  if (savePlot && !avgMeas) {
    ggsave(paste0("across_subj_mean_abs_bar_standarderror.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    ggsave(paste0("across_subj_mean_abs_bar_standarderror.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
    
  } else if (savePlot && avgMeas){
    ggsave(paste0("across_subj_mean_abs_bar_standarderror_avg.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    ggsave(paste0("across_subj_mean_abs_bar_standarderror_avg.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
    
  }
  
  
  # bar plot 
  figWidth = 8
  figHeight = 4
  goodVecBlock <- c("A/B 25","A/B 200","A/A 200")
  p14 <- ggplot(data =  dataListSummarizeBar, aes(x = blockType, y = meanAbs,fill=blockType)) +
    geom_bar(stat="identity")  +
    labs(x = expression(paste("Experimental Condition")),y=expression(paste("Absolute Difference from Baseline Peak-To-Peak (",mu,"V)")),fill="Experimental Condition",title = paste0("EP Difference from Baseline by Conditioning Protocol")) +
    scale_fill_brewer(palette="Dark2")
  print(p14)
  
  if (savePlot && !avgMeas) {
    ggsave(paste0("across_subj_mean_abs_bar.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    ggsave(paste0("across_subj_mean_abs_bar.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
    
  } else if (savePlot && avgMeas){
    ggsave(paste0("across_subj_mean_abs_bar_avg.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    ggsave(paste0("across_subj_mean_abs_bar_avg.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
    
  }
  
  
  # bar with dots
  figWidth = 8
  figHeight = 4
  goodVecBlock <- c("A/B 25","A/B 200","A/A 200")
  p13 <- ggplot() +
    geom_jitter(data =  dataListSummarize%>% filter(blockType %in% goodVecBlock),shape=16, position=position_jitter(0.2), aes(x = blockType, y = meanAbs,color=blockType,alpha = mapStimLevel)) +
    guides(fill=FALSE) +
    geom_bar(data=dataListSummarizeBar,stat="identity",aes(x = blockType, y = meanAbs,fill=blockType))  +
    labs(x = expression(paste("Experimental Condition")),y=expression(paste("Absolute Difference from Baseline Peak-To-Peak (",mu,"V)")),color="Experimental Condition",alpha="Ordered Stim Level",title = paste0("EP Difference from Baseline by Conditioning Protocol")) +
    scale_color_brewer(palette="Dark2") +
    scale_fill_brewer(palette="Dark2")
  print(p13)
  
  if (savePlot && !avgMeas) {
    ggsave(paste0("across_subj_mean_abs_bar_dots.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    ggsave(paste0("across_subj_mean_abs_bar_dots.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
    
  } else if (savePlot && avgMeas){
    ggsave(paste0("across_subj_mean_abs_bar_dots_avg.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    ggsave(paste0("across_subj_mean_abs_bar_dots_avg.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
    
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
    ggsave(paste0("diff_scale_each_subj_abs.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    ggsave(paste0("diff_scale_each_subj_abs.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
    
  } else if (savePlot && avgMeas){
    ggsave(paste0("diff_scale_each_subj_abs_avg.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    ggsave(paste0("diff_scale_each_subj_abs_avg.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
    
  }
  if (repeatedMeasures){
  fit.lmmPP = lmerTest::lmer(PPvec ~ mapStimLevel + disease  + chanInCond + overallBlockType*pre_post + baLabel + (1|subjectNum),data=dataList)
  emmeans(fit.lmmPP, list(pairwise ~ overallBlockType), adjust = "tukey")
  
  emm_s.t <- emmeans(fit.lmmPP, pairwise ~ overallBlockType| mapStimLevel)
  emm_s.t <- emmeans(fit.lmmPP, pairwise ~ mapStimLevel | overallBlockType)
  summary(glht(fit.lmmPP,linfct=mcp(overallBlockType="Tukey")))
  
  emm_pairwise <- emmeans(fit.lmmPP,~overallBlockType*pre_post,adjust="Tukey")
  contrast(emm_pairwise,interaction="pairwise")
  eff_size(emm_pairwise,sigma=sigma(fit.lmmPP),edf=df.residual(fit.lmmPP))
  marginal_means_plot <- emmip(fit.lmmPP,overallBlockType~pre_post)
  marginal_means_plot <- marginal_means_plot + aes(x = factor(pre_post, level=c('pre','post'))) + labs(x = expression(paste("Pre versus Post Conditioning")),y=expression(paste("Linear Prediction log(",mu,"V)")),color="Experimental Condition",title = paste0("Estimated Marginal Means by Conditioning and Status")) +
    scale_color_brewer(palette="Dark2") +
    scale_fill_brewer(palette="Dark2") 

  
  if (savePlot && !avgMeas) {
    ggsave(paste0("emmip.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    ggsave(paste0("emmip.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
    
  } else if (savePlot && avgMeas){
    ggsave(paste0("emmip_avg.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    ggsave(paste0("emmip_avg.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
    
  }
  
  }
  else {
  fit.lmmPP = lmerTest::lmer(PPvec ~ mapStimLevel + disease  + chanInCond + blockType + baLabel + (1|subjectNum),data=dataList)
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
  emmeans(fit.lmmPP, list(pairwise ~ baLabel), adjust = "tukey")
  

  anova(fit.lmmPP)
  tab_model(fit.lmmPP)

  
  if (savePlot && !avgMeas) {
    
    figHeight = 4
    figWidth = 8
    png("pairedPulse_resid_PP_trim_allSubjs.png",width=figWidth,height=figHeight,units="in",res=600)
    plot(fit.lmmPP)
    dev.off()
    
    setEPS()
    postscript("pairedPulse_resid_PP_trim_allSubjs.eps",width=figWidth,height=figHeight)
    plot(fit.lmmPP)
    dev.off()
    
    figHeight = 4
    figWidth = 8
    png("pairedPulse_qq_PP__trim_allSubjs.png",width=figWidth,height=figHeight,units="in",res=600)
    qqnorm(resid(fit.lmmPP))
    qqline(resid(fit.lmmPP))  #summary(fit.lmm2)
    dev.off()
    
    setEPS()
    postscript("pairedPulse_qq_PP_trim_allSubjs.eps",width=figWidth,height=figHeight)
    qqnorm(resid(fit.lmmPP))
    qqline(resid(fit.lmmPP))  #summary(fit.lmm2)
    dev.off()
    
  } else if (savePlot && avgMeas){
    
    figHeight = 4
    figWidth = 8
    png("pairedPulse_resid_PP_allSubjs_avg.png",width=figWidth,height=figHeight,units="in",res=600)
    plot(fit.lmmPP)
    dev.off()
    
    setEPS()
    postscript("pairedPulse_resid_PP_allSubjs_avg.eps",width=figWidth,height=figHeight)
    plot(fit.lmmPP)
    dev.off()
    
    figHeight = 4
    figWidth = 8
    png("pairedPulse_qq_PP_allSubjs_avg.png",width=figWidth,height=figHeight,units="in",res=600)
    qqnorm(resid(fit.lmmPP))
    qqline(resid(fit.lmmPP))  #summary(fit.lmm2)
    dev.off()
    
    setEPS()
    postscript("pairedPulse_qq_PP_allSubjs_avg.eps",width=figWidth,height=figHeight)
    qqnorm(resid(fit.lmmPP))
    qqline(resid(fit.lmmPP))  #summary(fit.lmm2)
    dev.off()
    
  }
  
  #fit.lmmdiff = lmerTest::lmer(absDiff ~ stimLevelVec + blockType + chanInCond + (1|sidVec),data=dataList)
  fit.lmmdiff = lmerTest::lmer(absDiff ~ mapStimLevel + disease  + chanInCond + blockType + (1|subjectNum),data=dataList)
  
  
  fit.lm = lm(absDiff ~ mapStimLevel + blockType + chanInCond ,data=dataList)
  
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
    png("pairedPulse_resid_diff_allSubjs.png",width=figWidth,height=figHeight,units="in",res=600)
    plot(fit.lmmdiff)
    dev.off()
    
    setEPS()
    postscript("pairedPulse_resid_diff_allSubjs.eps",width=figWidth,height=figHeight)
    plot(fit.lmmdiff)
    dev.off()
    
    figHeight = 4
    figWidth = 8
    png("pairedPulse_qq_diff_allSubjs.png",width=figWidth,height=figHeight,units="in",res=600)
    qqnorm(resid(fit.lmmdiff))
    qqline(resid(fit.lmmdiff))  #summary(fit.lmm2)
    dev.off()
    
    setEPS()
    postscript("pairedPulse_qq_diff_allSubjs.eps",width=figWidth,height=figHeight)
    qqnorm(resid(fit.lmmdiff))
    qqline(resid(fit.lmmdiff))  #summary(fit.lmm2)
    dev.off()
    
  } else if (savePlot && avgMeas){
    
    figHeight = 4
    figWidth = 8
    png("pairedPulse_qq_resid_diff_allSubjs_avg.png",width=figWidth,height=figHeight,units="in",res=600)
    plot(fit.lmmdiff)
    dev.off()
    
    setEPS()
    postscript("pairedPulse_resid_diff_allSubjs_avg.eps",width=figWidth,height=figHeight)
    plot(fit.lmmdiff)
    dev.off()
    
    figHeight = 4
    figWidth = 8
    png("pairedPulse_qq_diff_allSubjs_avg.png",width=figWidth,height=figHeight,units="in",res=600)
    qqnorm(resid(fit.lmmdiff))
    qqline(resid(fit.lmmdiff))  #summary(fit.lmm2)
    dev.off()
    
    setEPS()
    postscript("pairedPulse_qq_diff_allSubjs_avg.eps",width=figWidth,height=figHeight)
    qqnorm(resid(fit.lmmdiff))
    qqline(resid(fit.lmmdiff))  #summary(fit.lmm2)
    dev.off()
    
  }
  
  dataSubset <- unique(dataList %>% select(effectSize,subjectNum,meanPP,mapStimLevel,disease,chanInCond,blockType,baLabel,aalLabel) %>% filter(blockType != 'baseline'))
  if(!log_data){
  fit.effectSize = lmerTest::lmer(effectSize ~ log(meanPP) + mapStimLevel + chanInCond + disease + blockType + baLabel + (1|subjectNum),data=dataSubset)
  }else if(log_data){
    fit.effectSize = lmerTest::lmer(effectSize ~ meanPP + mapStimLevel + chanInCond + disease + blockType + baLabel + (1|subjectNum),data=dataSubset)
  }
  
  # if(!log_data){
  #   fit.effectSize = lmerTest::lmer(effectSize ~ log(meanPP) + chanInCond + disease + blockType + (1|subjectNum),data=dataSubset)
  # }else if(log_data){
  #   fit.effectSize = lmerTest::lmer(effectSize ~ meanPP + chanInCond + disease + blockType + (1|subjectNum),data=dataSubset)
  # }
  # 
  emmeans(fit.effectSize, list(pairwise ~ blockType), adjust = "tukey")
  
  
  
  #fit.effectSize = lm(effectSize ~ meanPP + mapStimLevel + chanInCond + disease + blockType ,data=dataSubset)
  
  plot(fit.effectSize)
  
  qqnorm(resid(fit.effectSize))
  qqline(resid(fit.effectSize))  #summary(fit.lmm2)
  
  figWidth = 8
  figHeight = 4
  goodVecBlock <- c("A/B 25","A/B 200","A/A 200")
  pEffect <- ggplot(data =  dataSubset%>% filter(blockType %in% goodVecBlock), aes(x = blockType, y = effectSize,color=blockType)) +
    geom_boxplot(notch=TRUE,outlier.shape=NA)  + geom_jitter(shape=16, position=position_jitter(0.2),aes(alpha = mapStimLevel)) +
    labs(x = expression(paste("Experimental Condition")),y=expression(paste("Effect Size")),color="Experimental Condition",alpha="Ordered Stim Level",title = paste0("Effect Size by Conditioning Protocol")) +
    scale_color_brewer(palette="Dark2") + scale_alpha_discrete(range=c(0.5,1))
  print(pEffect)
  if (savePlot && !avgMeas) {
    ggsave(paste0("across_subj_effect_trim_all_subjs.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    ggsave(paste0("across_subj_effect_trim_all_subjs.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
    
  } else if (savePlot && avgMeas){
    ggsave(paste0("across_subj_effect_trim_avg.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    ggsave(paste0("across_subj_effect_trim_avg.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
    
  }
  
  if (savePlot && !avgMeas) {
    
    figHeight = 4
    figWidth = 8
    png("pairedPulse_resid_effectSize_allSubjs.png",width=figWidth,height=figHeight,units="in",res=600)
    plot(fit.effectSize)
    dev.off()
    
    setEPS()
    postscript("pairedPulse_resid_effectSize_allSubjs.eps",width=figWidth,height=figHeight)
    plot(fit.effectSize)
    dev.off()
    
    figHeight = 4
    figWidth = 8
    png("pairedPulse_qq_effectSize_allSubjs.png",width=figWidth,height=figHeight,units="in",res=600)
    qqnorm(resid(fit.effectSize))
    qqline(resid(fit.effectSize))  #summary(fit.lmm2)
    dev.off()
    
    setEPS()
    postscript("pairedPulse_qq_effectSize_allSubjs.eps",width=figWidth,height=figHeight)
    qqnorm(resid(fit.effectSize))
    qqline(resid(fit.effectSize))  #summary(fit.lmm2)
    dev.off()
    
  } else if (savePlot && avgMeas){
    
    figHeight = 4
    figWidth = 8
    png("pairedPulse_resid_effectSize_allSubjs_avg.png",width=figWidth,height=figHeight,units="in",res=600)
    plot(fit.effectSize)
    dev.off()
    
    setEPS()
    postscript("pairedPulse_resid_effectSize_allSubjs_avg.eps",width=figWidth,height=figHeight)
    plot(fit.lmmPP)
    dev.off()
    
    figHeight = 4
    figWidth = 8
    png("pairedPulse_qq_effectSize_allSubjs_avg.png",width=figWidth,height=figHeight,units="in",res=600)
    qqnorm(resid(fit.lmmPP))
    qqline(resid(fit.lmmPP))  #summary(fit.lmm2)
    dev.off()
    
    setEPS()
    postscript("pairedPulse_qq_effectSize_allSubjs_avg.eps",width=figWidth,height=figHeight)
    qqnorm(resid(fit.lmmPP))
    qqline(resid(fit.lmmPP))  #summary(fit.lmm2)
    dev.off()
    
  }
  
}

