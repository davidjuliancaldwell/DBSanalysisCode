
library('here')
library('nlme')
library('ggplot2')
library('drc')
library('minpack.lm')
rootDir = here()
dataDir = here("DBS_CCEP_PairedPulse","R_data")
codeDir = here("DBS_CCEP_PairedPulse","R_code")

#sidVec <- c('46c2a','08b13','8e907')
sidVec <- c('46c2a')

for (sid in sidVec){
  source(here("DBS_CCEP_PairedPulse","R_config_files",paste0("subj_",sid,".R")))
  
  dataPP <- read.table(here("DBS_CCEP_PairedPulse","R_data",paste0(sid,'_PairedPulseData.csv')),header=TRUE,sep = ",",stringsAsFactors=F, colClasses=c("stimLevelVec"="numeric","sidVec"="character"))
  print(sid)
  for (chanInt in chanIntVec){
    dataInt <- na.exclude(subset(dataPP,chanVec == chanInt & blockVec %in% blockIntPlot))
    p <- ggplot(dataInt, aes(x=stimLevelVec, y=1e6*PPvec,colour=stimLevelVec)) + 
      geom_jitter(width=35) +  geom_smooth(method=lm) + facet_wrap(~blockVec, labeller = as_labeller(blockNames))+
      labs(x = expression(paste("Stimulation current ",mu,"A")),y=expression(paste("Peak to Peak magnitude ",mu,"V"), fill="stimulus level"),title = paste0("Subject ", subjectNum, " ID ", sid," DBS Paired Pulse EP Measurements")) + guides(colour=guide_legend("stimulation level")) 
    print(p)
    
    
    dataToFit <- na.exclude(dataPP[dataPP$blockVec %in% blockIntLM & dataPP$chanVec==chanInt,])
    fit.glm    = glm(PPvec ~ as.factor(blockVec) + stimLevelVec,data=dataToFit)
    summary(fit.glm)
    plot(fit.glm)
    
    fitnlme0 <- nlme(1e6*PPvec ~ SSlogis(stimLevelVec, Asym, xmid, scal),
                     fixed=list(Asym + xmid + scal ~ 1),
                    random = Asym ~ 1 | blockVec,
                     start=list(fixed=c(Asym=100,xmid=1500,scal=500)),
                     data=dataToFit)
    
    
    fitnlme1 <- nlme(1e6*PPvec ~ SSlogis(stimLevelVec, Asym, xmid, scal),
                     fixed=list(Asym ~ blockVec, xmid + scal ~ 1),
                   random = Asym ~ 1 | blockVec,
                     start=list(fixed=c(Asym=c(100,140),xmid=1500,scal=500)),
                     data=dataToFit)
    
    # Likelihood ratio test
    anova(fitnlme0, fitnlme1)  
  }
}

#http://rstudio-pubs-static.s3.amazonaws.com/28730_850cac53898b45da8050f7f622d48927.html