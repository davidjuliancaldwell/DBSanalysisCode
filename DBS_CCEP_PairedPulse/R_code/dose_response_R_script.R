
library('here')
library('nlme')
library('ggplot2')
library('drc')
library('minpack.lm')
library('lmtest')
rootDir = here()
dataDir = here("DBS_CCEP_PairedPulse","R_data")
codeDir = here("DBS_CCEP_PairedPulse","R_code")

#theme_set(theme_grey(base_size = 24)) 

sidVec <- c('46c2a','9f852','08b13','8e907')

savePlot = FALSE

for (sid in sidVec){
  source(here("DBS_CCEP_PairedPulse","R_config_files",paste0("subj_",sid,".R")))
  
  dataPP <- read.table(here("DBS_CCEP_PairedPulse","R_data",paste0(sid,'_PairedPulseData_avg.csv')),header=TRUE,sep = ",",stringsAsFactors=F, colClasses=c("stimLevelVec"="numeric","sidVec"="character"))
  # multiply by 1e6
  dataPP$PPvec = dataPP$PPvec*1e6
  # change to factor 
  dataPP$blockVec = as.factor(dataPP$blockVec)
  print(sid)
  for (chanInt in chanIntVec){
    
    dataInt <- na.exclude(subset(dataPP,(chanVec == chanInt) & (blockVec %in% blockIntPlot)))
    
    p <- ggplot(dataInt, aes(x=stimLevelVec, y=PPvec,colour=stimLevelVec)) + 
      geom_jitter(width=35) +  geom_smooth(method=lm) + facet_wrap(~blockVec, labeller = as_labeller(blockNames))+
      labs(x = expression(paste("Stimulation current ",mu,"A")),y=expression(paste("Peak to Peak magnitude ",mu,"V"), fill="stimulus level"),title = paste0("Subject ", subjectNum, " ID ", sid," DBS Paired Pulse EP Measurements")) +
      guides(colour=guide_legend("stimulation level"))
    print(p)
    
    if (savePlot) {
    ggsave(paste0("subj_", subjectNum, "_ID_", sid,"_scatter_lm.png"), units="in", width=6, height=4, dpi=600)
    }
    
    dataToFit <- na.exclude(dataPP[dataPP$blockVec %in% blockIntLM & dataPP$chanVec==chanInt,])
    fit.glm    = glm(PPvec ~ blockVec + stimLevelVec,data=dataToFit)
    summary(fit.glm)
   # plot(fit.glm)
    
  #  fitnlme0 <- nlme(PPvec ~ SSlogis(stimLevelVec, Asym, xmid, scal),
  #                   fixed=list(Asym + xmid + scal ~ 1),
   #                  random = Asym ~ 1 | blockVec,
  #                   start=list(fixed=c(Asym=100,xmid=1500,scal=500)),
   #                  data=dataToFit)
    
    # get starting values 
    fit.startnlme0 <- nlsLM(PPvec ~ SSlogis(stimLevelVec, Asym, xmid, scal),data=dataToFit)
    startVals = summary(fit.startnlme0)$coefficients
    fit.nlme0 <- gnls(PPvec ~ SSlogis(stimLevelVec, Asym, xmid, scal),data=dataToFit,
                    start=list(Asym=startVals[1,1]+400,xmid=startVals[2,1],scal=startVals[3,1]))
    
    # fit.nlme0 <- gnls(PPvec ~ SSlogis(stimLevelVec, Asym, xmid, scal),data=dataToFit,
    #                   start=list(Asym=1116,xmid=3921,scal=1062))
    # 
    # 
    
 
    
#     existData = data.frame(PPvec = predict(fit.nlme0,newData=dataToFit),stimLevelVec = dataToFit$stimLevelVec)
# 
#     p2 <- ggplot(dataToFit, aes(x=stimLevelVec, y=PPvec,colour=stimLevelVec)) +
#       geom_jitter(width=35) +
#       facet_wrap(~blockVec, labeller = as_labeller(blockNames))+
#       labs(x = expression(paste("Stimulation current ",mu,"A")),y=expression(paste("Peak to Peak magnitude ",mu,"V"), fill="stimulus level"),title = paste0("Subject ", subjectNum, " ID ", sid," DBS Paired Pulse EP Measurements")) +
#       guides(colour=guide_legend("stimulation level")) +
#       geom_line(data=existData,aes(x=stimLevelVec, y=PPvec))
# print(p2)


# p2 <- ggplot(dataInt, aes(x=stimLevelVec, y=1e6*PPvec,colour=stimLevelVec)) + 
#   geom_jitter(width=35) +  geom_smooth(method="nls",formula=y~Asym/(1+exp((xmid-x)/scal)),method.args=list(start=c(Asym=443,xmid=3990,scal=900)),se=FALSE) +
#   facet_wrap(~blockVec, labeller = as_labeller(blockNames))+
#   labs(x = expression(paste("Stimulation current ",mu,"A")),y=expression(paste("Peak to Peak magnitude ",mu,"V"), fill="stimulus level"),title = paste0("Subject ", subjectNum, " ID ", sid," DBS Paired Pulse EP Measurements")) + 
#   guides(colour=guide_legend("stimulation level"))+ geom_line(data=predData,aes(x=stimLevelVec,y=PPvec))
# print(p2)
# 

#ggplot(dataInt,aes(x=factor(stimLevelVec),y=PPvec)) + geom_violin() + geom_line(data=predData,inherit.aes = FALSE,aes(y=PPvec,x=factor(stimLevelVecNew$stimLevelVec),group=1))

#Create a range of doses:
stimLevelVecNew <- data.frame(stimLevelVec = seq(min(dataToFit$stimLevelVec), max(dataToFit$stimLevelVec), length.out = 200))
#Create a new data frame for ggplot using predict and your range of new 
#doses:
predData=data.frame(PPvec=predict(fit.nlme0,newdata=stimLevelVecNew),stimLevelVec=stimLevelVecNew)

p3 <- ggplot(dataToFit, aes(x=stimLevelVec, y=PPvec,colour=stimLevelVec)) + 
  geom_jitter(width=35) +
  facet_wrap(~blockVec, labeller = as_labeller(blockNames))+
  labs(x = expression(paste("Stimulation current ",mu,"A")),y=expression(paste("Peak to Peak magnitude ",mu,"V"), fill="stimulus level"),title = paste0("Subject ", subjectNum, " ID ", sid," DBS Paired Pulse EP Measurements")) + 
  guides(colour=guide_legend("stimulation level"))+ geom_line(data=predData,aes(x=stimLevelVec,y=PPvec),size=2,colour='black')
print(p3)

if (savePlot) {
ggsave(paste0("subj_", subjectNum, "_ID_", sid,"_scatter_sameModel.png"), units="in", width=6, height=4, dpi=600)
}

fit.nlme1 <- nlme(PPvec ~ SSlogis(stimLevelVec, Asym, xmid, scal),
                 fixed=list(Asym ~ blockVec, xmid + scal ~ 1),
                 random = scal ~ 1|blockVec,
                 start=list(fixed=c(Asym=rep(startVals[1,1]+200,length(blockIntLM)),xmid=startVals[2,1],scal=startVals[3,1])),
                 data=dataToFit)

# vs.   fixed=list(Asym ~ blockVec, xmid + scal ~ 1),


#Create a range of doses:
predDataGroups <-expand.grid(stimLevelVec = seq(min(dataToFit$stimLevelVec), max(dataToFit$stimLevelVec), length.out = 100),blockVec = unique(dataToFit$blockVec))
#Create a new data frame for ggplot using predict and your range of new 
#doses:
predDataGroups$PPvec=predict(fit.nlme1,newdata=predDataGroups)

p4 <- ggplot(dataToFit, aes(x=stimLevelVec, y=PPvec,colour=stimLevelVec)) + 
  geom_jitter(width=35) +
  facet_wrap(~blockVec, labeller = as_labeller(blockNames))+
  labs(x = expression(paste("Stimulation current ",mu,"A")),y=expression(paste("Peak to Peak magnitude ",mu,"V"), fill="stimulus level"),title = paste0("Subject ", subjectNum, " ID ", sid," DBS Paired Pulse EP Measurements")) + 
  guides(colour=guide_legend("stimulation level"))+ 
  geom_line(data=predDataGroups,aes(x=stimLevelVec,y=PPvec),size=2,colour='black')

print(p4)

if (savePlot) {
ggsave(paste0("subj_", subjectNum, "_ID_", sid,"_scatter_diffModel.png"), units="in", width=6, height=4, dpi=600)
}

# xtra ss test
anova(fit.nlme0, fit.nlme1)  
# Likelihood ratio test

lrtest(fit.nlme0, fit.nlme1,fit.glm)  

# aic
AIC(fit.nlme0, fit.nlme1,fit.glm)  

  }
}

#http://rstudio-pubs-static.s3.amazonaws.com/28730_850cac53898b45da8050f7f622d48927.html