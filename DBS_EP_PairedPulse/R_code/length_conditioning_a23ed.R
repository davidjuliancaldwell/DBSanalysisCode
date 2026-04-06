
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
library('dplyr')
library('DescTools')
library('effectsize')

rootDir = here()
dataDir = here("DBS_EP_PairedPulse","R_data")
codeDir = here("DBS_EP_PairedPulse","R_code")
outputDir = here("DBS_EP_PairedPulse","R_output")
dir.create(outputDir, showWarnings = FALSE)

sidVec = c("a23ed")


repeatedMeasures = TRUE # if true, does repeated measures analysis, if false, does more of ANCOVA style analysis
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

# hardcode subject #
subjectNum = 11

oldwd <- getwd()
setwd(outputDir)

for (avgMeas in avgMeasVec) {
  for (sid in sidVec){
    source(here("DBS_EP_PairedPulse","R_config_files",paste0("subj_",sid,".R")))
    
    # just modifiy a23ed
    whichCompareVec = list(c(2,3),c(5,6))
    blockType = c('baseline','baseline','A/B 200 ms 5 minutes','baseline','baseline',
                  'A/B 200 ms 15 minutes','baseline','A/A 200','baseline')
    
    dataPP <- read.table(here("DBS_EP_PairedPulse","R_data",paste0(sid,'_PairedPulseData_avg_5.csv')),header=TRUE,sep = ",",stringsAsFactors=F, colClasses=c("stimLevelVec"="numeric","sidVec"="character"))

    # PPvec is in volts, convert to microvolts
    dataPP$PPvec = dataPP$PPvec*1e6
    dataPP$stimLevelVec = dataPP$stimLevelVec/1e3
    #dataPP <- subset(dataPP, PPvec<1000)
    #dataPP <- subset(dataPP, PPvec>30)
    
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
      dataInt <- subset(dataPP,(chanVec == chanInt) & (blockVec %in% blockIntPlot))
      
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
      
      #uniqueBlockLevel = unique(dataInt$blockVec)
      #blockTypeTrim = blockType[uniqueBlockLevel]
      
      #dataInt$blockType <- mapvalues(dataInt$blockVec,
      #                               from=uniqueBlockLevel,
      #                               to=blockTypeTrim)
      
      
      dataInt$blockType <- mapvalues(dataInt$blockVec,
                                     from=blockIntPlot,
                                     to=blockType)
      
      dataInt$mapStimLevel = as.ordered(dataInt$mapStimLevel)
      dataInt$blockType = as.factor(dataInt$blockType)
      
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
        
        dataIntCompare <- dataIntCompare %>% 
          group_by(mapStimLevel) %>% 
          mutate(absDiff = PPvec - mean(PPvec[blockVec==comparison[1]]),
                 percentDiff = 100*(PPvec - mean(PPvec[blockVec==comparison[1]]))/mean(PPvec[blockVec==comparison[1]]) ) 
        
        dataIntCompare = as_data_frame(dataIntCompare)
        dataIntCompare$index = index
        
        
        name <- unique(dataIntCompare %>% filter(blockType!='baseline')%>%select(blockType))
        name <- as.character(name$blockType)
        
        dataIntCompare$overallBlockType <- as.factor(name)
        dataIntCompare <- dataIntCompare %>% mutate(pre_post = case_when(blockType == name ~ 'post',
                                                                         blockType=='baseline' ~'pre'
        ))
        dataList[[index]] = dataIntCompare
        blockList[[index]] = comparison
        
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
        
        index = index + 1
        
      }
      
      p <- ggplot(dataInt, aes(x=as.numeric(mapStimLevel), y=PPvec,color=stimLevelVec)) +
        geom_point(position=position_jitterdodge(dodge.width=0.250)) +  geom_smooth(method=lm) + facet_wrap(~blockVec, labeller = as_labeller(blockNames))+
        labs(x = expression(paste("Stimulation Current (mA)")),y=expression(paste("Peak to Peak Magnitude (",mu,"V)")),title = paste0("Subject ", subjectNum," DBS Paired Pulse EP Measurements")) +
        guides(colour=guide_colorbar("Stimulation Level"))
      print(p)
      
      if (savePlot && !avgMeas) {
        ggsave(paste0("subj_", subjectNum, "_ID_", sid,"_chan_",chanInt,"_compare_length_scatter_lm.png"), units="in", width=figWidth, height=figHeight, dpi=600)
        ggsave(paste0("subj_", subjectNum, "_ID_", sid,"_chan_",chanInt,"_compare_length_scatter_lm.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
        
      }
      else if (savePlot && avgMeas){
        ggsave(paste0("subj_", subjectNum, "_ID_", sid,"_chan_",chanInt,"_compare_length_scatter_lm_avg.png"), units="in", width=figWidth, height=figHeight, dpi=600)
        ggsave(paste0("subj_", subjectNum, "_ID_", sid,"_chan_",chanInt,"_compare_length_scatter_lm_avg.eps"), units="in", width=figWidth, height=figHeight, device=cairo_ps, fallback_resolution=600)
        
      }
      
      dataToFit <- na.exclude(dataPP[dataPP$blockVec %in% blockIntLM & dataPP$chanVec==chanInt,])
      fit.glm    = glm(PPvec ~ blockVec + stimLevelVec,data=dataToFit)
      summary(fit.glm)
      summary(glht(fit.glm,linfct=mcp(blockVec="Tukey")))
      
    }
    
    
    # now do comparison at highest stim level relative to baseline
    dataList <- do.call(rbind,dataList)
    #blockList <- do.call(rbind,blockList)

    # Aggregate to cell medians: one row per (block, stim_level)
    # to avoid trial-level pseudoreplication
    dataListAgg <- dataList %>%
      group_by(chanVec, blockVec, blockType, mapStimLevel, stimLevelVec,
               overallBlockType, pre_post, sidVec, index) %>%
      summarise(PPvec = median(PPvec), absDiff = median(absDiff),
                percentDiff = median(percentDiff), n_trials = n(), .groups = "drop")

    cat("Aggregated:", nrow(dataListAgg), "cell medians from", nrow(dataList), "trials\n")

    #plot
    grouped <- group_by(dataList, sidVec, chanVec, blockType,mapStimLevel)
    dataListSummarize <- summarise(grouped,meanPerc = mean(percentDiff),sdPerc = sd(percentDiff),
                                   meanAbs=mean(absDiff), sdDiff=sd(percentDiff),meanPP = mean(PPvec),sdPP = sd(PPvec))
    
    
    p3 <- ggplot(dataList, aes(x=as.numeric(mapStimLevel), y=percentDiff,color=blockType)) +
      geom_point(position=position_jitterdodge(dodge.width=0.75)) +
      labs(x = expression(paste("Stimulation level")),y=expression(paste("Percent Difference in EP Magnitude from Baseline")),
           title = paste0("Changes in EP magnitude"),color="Experimental Condition")
    print(p3)
    
    if (savePlot && !avgMeas) {
      ggsave(paste0("across_a23ed_percent.png"), units="in", width=figWidth, height=figHeight, dpi=600)
      ggsave(paste0("across_a23ed_percent.eps"), units="in", width=figWidth, height=figHeight, device=cairo_ps, fallback_resolution=600)
      
    } else if (savePlot && avgMeas){
      ggsave(paste0("across_a23ed_percent_avg.png"), units="in", width=figWidth, height=figHeight, dpi=600)
      ggsave(paste0("across_a23ed_percent_avg.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
      
    }
    
    p4 <- ggplot(dataList, aes(x=as.numeric(mapStimLevel), y=absDiff,color=blockType)) +
      geom_point(position=position_jitterdodge(dodge.width=0.75)) +geom_smooth(method=lm) +
      labs(x = expression(paste("Stimulation Level")),y=expression(paste("Absolute Difference in EP Magnitude from Baseline (",mu,"V)")), color="Experimental Condition",title = paste0("Changes in EP magnitude"))
    print(p4)
    
    if (savePlot && !avgMeas) {
      ggsave(paste0("across_a23ed_abs.png"), units="in", width=figWidth, height=figHeight, dpi=600)
      ggsave(paste0("across_a23ed_abs.eps"), units="in", width=figWidth, height=figHeight, device=cairo_ps, fallback_resolution=600)
      
    } else if (savePlot && avgMeas){
      ggsave(paste0("across_a23ed_abs_avg.png"), units="in", width=figWidth, height=figHeight, dpi=600)
      ggsave(paste0("across_a23ed_abs_avg.eps"), units="in", width=figWidth, height=figHeight, device=cairo_ps, fallback_resolution=600)
      
    }
    
    p5 <- ggplot(dataListSummarize, aes(x=as.numeric(mapStimLevel), y=meanPerc,color=blockType)) +
      geom_point(position=position_jitterdodge(dodge.width=0.5)) +geom_smooth(method=lm) +
      labs(x = expression(paste("Stimulation Level")),y=expression(paste("Mean Percent Difference in EP Magnitude")), color="Experimental Condition",
           title = paste0("Changes in EP Magnitude"))
    print(p5)
    
    if (savePlot && !avgMeas) {
      ggsave(paste0("across_a23ed_mean_perc.png"), units="in", width=figWidth, height=figHeight, dpi=600)
      ggsave(paste0("across_a23ed_mean_perc.eps"), units="in", width=figWidth, height=figHeight, device=cairo_ps, fallback_resolution=600)
      
    } else if (savePlot && avgMeas){
      ggsave(paste0("across_a23ed_mean_perc_avg.png"), units="in", width=figWidth, height=figHeight, dpi=600)
      ggsave(paste0("across_a23ed_mean_perc_avg.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
      
    }
    
    p6 <- ggplot(dataListSummarize, aes(x=as.numeric(mapStimLevel), y=meanAbs,color=blockType)) +
      geom_point(position=position_jitterdodge(dodge.width=0.5)) +geom_smooth(method=lm) +
      labs(x = expression(paste("Stimulation Level")),y=expression(paste("Mean Absolute Difference in EP Magnitude from Baseline (",mu,"V)")),
           color="Experimental Condition",title = paste0("Changes in EP Magnitude"))
    print(p6)
    
    if (savePlot && !avgMeas) {
      ggsave(paste0("across_a23ed_mean_abs.png"), units="in", width=figWidth, height=figHeight, dpi=600)
      ggsave(paste0("across_a23ed_mean_abs.eps"), units="in", width=figWidth, height=figHeight, device=cairo_ps, fallback_resolution=600)
      
    } else if (savePlot && avgMeas){
      ggsave(paste0("across_a23ed_mean_abs_avg.png"), units="in", width=figWidth, height=figHeight, dpi=600)
      ggsave(paste0("across_a23ed_mean_abs_avg.eps"), units="in", width=figWidth, height=figHeight, device=cairo_ps, fallback_resolution=600)
      
    }
    
    p7 <- ggplot(dataList, aes(x=as.numeric(mapStimLevel), y=PPvec,color=blockType)) +
      geom_point(position=position_jitterdodge(dodge.width=0.75)) +geom_smooth(method=lm) +
      labs(x = expression(paste("Stimulation Level")),y=expression(paste("Peak to Peak voltage (",mu,"V)")),
           color="Experimental Condition",title = paste0("EP Magnitude by Length of Conditioning"))
    print(p7)
    # if (savePlot && !avgMeas) {
    #   ggsave(paste0("across_a23ed_PP.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    #   ggsave(paste0("across_a23ed_PP.eps"), units="in", width=figWidth, height=figHeight, dpi=600,device=cairo_ps, fallback_resolution=600)
    #   
    # } else if (savePlot && avgMeas){
    #   ggsave(paste0("across_a23ed_PP_avg.png"), units="in", width=figWidth, height=figHeight, dpi=600)
    #   ggsave(paste0("across_a23ed_PP_avg.eps"), units="in", width=figWidth, height=figHeight, dpi=600,device=cairo_ps, fallback_resolution=600)
    #   
    #   
    # }
    # 
    # 
    
    
    p7 <- ggplot(dataList, aes(x=as.numeric(mapStimLevel), y=PPvec,color=blockType)) +
      geom_point(position=position_jitterdodge(dodge.width=0.75)) +geom_smooth(method=lm) +
      labs(x = expression(paste("Stimulation Level")),y=expression(paste("Peak to Peak voltage (",mu,"V)")),
           color="Experimental Condition",title = paste0("EP Magnitude by Length of Conditioning")) + scale_color_manual(values=c("#F8766D","#bedd73","#7CAE00" ))
    
    print(p7)
    if (savePlot && !avgMeas) {
      ggsave(paste0("across_a23ed_PP.png"), units="in", width=figWidth, height=figHeight, dpi=600)
      ggsave(paste0("across_a23ed_PP.eps"), units="in", width=figWidth, height=figHeight, dpi=600,device=cairo_ps, fallback_resolution=600)
      
    } else if (savePlot && avgMeas){
      ggsave(paste0("across_a23ed_PP_avg.png"), units="in", width=figWidth, height=figHeight, dpi=600)
      ggsave(paste0("across_a23ed_PP_avg.eps"), units="in", width=figWidth, height=figHeight, dpi=600,device=cairo_ps, fallback_resolution=600)
      
    }
    
    
    p8 <- ggplot(dataListSummarize, aes(x=as.numeric(mapStimLevel), y=meanPP,colour=blockType,fill=blockType)) +
      geom_point(position=position_jitterdodge(dodge.width=0.5)) +geom_smooth(method=lm) +
      labs(x = expression(paste("Stimulation Level")),y=expression(paste("Mean Peak to Peak Voltage (",mu,"V)")), color="Experimental Condition",
           title = paste0("EP Magnitude by Length of Conditioning"))
    print(p8)
    
    if (savePlot && !avgMeas) {
      ggsave(paste0("across_a23ed_mean_PP.png"), units="in", width=figWidth, height=figHeight, dpi=600)
      ggsave(paste0("across_a23ed_mean_PP.eps"), units="in", width=figWidth, height=figHeight, dpi=600,device=cairo_ps, fallback_resolution=600)
      
    } else if (savePlot && avgMeas){
      ggsave(paste0("across_a23ed_mean_PP_avg.png"), units="in", width=figWidth, height=figHeight, dpi=600)
      ggsave(paste0("across_a23ed_mean_PP_avg.eps"), units="in", width=figWidth, height=figHeight, dpi=600,device=cairo_ps, fallback_resolution=600)
      
    }  
    
    if (repeatedMeasures){
      # --- Model 1: Full interaction (overallBlockType * pre_post) ---
      # Tests whether 5-min vs 15-min conditioning differ in pre->post change.
      # NOTE: singular fit expected -- 4 fixed-effect combos saturate 4 blocks.
      # Included for completeness; effect size plots are primary for this comparison.
      fit.lmmPP.full = lmerTest::lmer(PPvec ~ mapStimLevel + overallBlockType*pre_post + (1|blockVec), data=dataList)
      cat("\n=== Model 1: Full interaction (may be singular) ===\n")
      print(summary(fit.lmmPP.full))

      # --- Model 2: pre_post only (no overallBlockType) ---
      # Tests whether conditioning changes EP overall (pooling across durations).
      # Block RE is estimable here (2 blocks per pre_post level).
      fit.lmmPP.prepost = lmerTest::lmer(PPvec ~ mapStimLevel + pre_post + (1|blockVec), data=dataList)
      cat("\n=== Model 2: pre_post only ===\n")
      print(summary(fit.lmmPP.prepost))

      cat("\nAIC: full =", AIC(fit.lmmPP.full), " pre_post =", AIC(fit.lmmPP.prepost), "\n")

      # Use pre_post model as primary (non-singular, proper block RE)
      fit.lmmPP = fit.lmmPP.prepost

      emmeans(fit.lmmPP, list(pairwise ~ pre_post), adjust = "tukey")

      emm_pairwise <- emmeans(fit.lmmPP, ~pre_post, adjust="tukey")
      contrast(emm_pairwise, method="pairwise")
      eff_size(emm_pairwise, sigma=sigma(fit.lmmPP), edf=df.residual(fit.lmmPP))

      emm_s.t <- emmeans(fit.lmmPP, pairwise ~ mapStimLevel)

      # Also show full model emmeans for the 5-min vs 15-min comparison
      emm_full <- emmeans(fit.lmmPP.full, ~overallBlockType*pre_post, adjust="tukey")
      contrast(emm_full, interaction="pairwise")

      marginal_means_plot <- emmip(fit.lmmPP.full, overallBlockType~pre_post)
      marginal_means_plot + labs(x = "Pre versus Post Conditioning",
        y = expression(paste("Linear Prediction (",mu,"V)")),
        color = "Experimental Condition",
        title = "Estimated Marginal Means by Conditioning Length and Pre versus Post")

      
      if (savePlot && !avgMeas) {
        ggsave(paste0("emmip_a23ed.png"), units="in", width=figWidth, height=figHeight, dpi=600)
        ggsave(paste0("emmip_a23ed.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
        
      } else if (savePlot && avgMeas){
        ggsave(paste0("emmip_avg_a23ed.png"), units="in", width=figWidth, height=figHeight, dpi=600)
        ggsave(paste0("emmip_avg_a23ed.eps"), units="in", width=figWidth, height=figHeight,device=cairo_ps, fallback_resolution=600)
        
      }

      }
    else{

    fit.lmmPP = lmerTest::lmer(PPvec ~ mapStimLevel + blockType + (1|blockVec), data=dataList)
    emmeans(fit.lmmPP, list(pairwise ~ blockType), adjust = "tukey")
    
    emm_s.t <- emmeans(fit.lmmPP, pairwise ~ blockType | mapStimLevel)
    
    emm_s.t <- emmeans(fit.lmmPP, pairwise ~ mapStimLevel | blockType)
    summary(glht(fit.lmmPP,linfct=mcp(blockType="Tukey")))
    
    }
    summary(fit.lmmPP)
    # plot(fit.lm)
    
    
    summary(glht(fit.lmmPP,linfct=mcp(mapStimLevel="Tukey")))
    
    
    emmeans(fit.lmmPP, list(pairwise ~ mapStimLevel), adjust = "tukey")
  
    
    anova(fit.lmmPP)
    tab_model(fit.lmmPP)
    tab_a23ed <- tab_model(fit.lmmPP.full, fit.lmmPP.prepost,
      show.re.var=TRUE, show.icc=TRUE, show.obs=TRUE,
      dv.labels=c("Full (overallBlockType*pre_post)","Reduced (pre_post only)"))
    writeLines(as.character(tab_a23ed$page.complete), paste0(outputDir, "/tab_model_a23ed_comparison.html"))
    
    fit.lmmdiff = lmerTest::lmer(absDiff ~ mapStimLevel + blockType + (1|blockVec), data=dataList)

    fit.lm = lm(absDiff ~ mapStimLevel + blockType, data=dataList)
    tab_model(fit.lm)
    
    summary(fit.lmmdiff)
    # plot(fit.lm)
    plot(fit.lmmPP)
    qqnorm(resid(fit.lmmPP))
    qqline(resid(fit.lmmPP))
    summary(glht(fit.lmmdiff,linfct=mcp(blockType="Tukey")))
    summary(glht(fit.lmmdiff,linfct=mcp(mapStimLevel="Tukey")))
    
    emmeans(fit.lmmdiff, list(pairwise ~ blockType), adjust = "tukey")
    emmeans(fit.lmmdiff, list(pairwise ~ mapStimLevel), adjust = "tukey")
    
    emm_s.t <- emmeans(fit.lmmdiff, pairwise ~ blockType | mapStimLevel)
    emm_s.t <- emmeans(fit.lmmdiff, pairwise ~ mapStimLevel | blockType)
    
    anova(fit.lmmdiff)
    tab_model(fit.lmmdiff)

    #########
    # Effect size analysis (appropriate for single-subject design)
    # Cohen's d per (stim_level, comparison): post-conditioning vs baseline
    # comparison[1] = baseline block, comparison[2] = post-conditioning block
    # NOTE: log_data=FALSE for this script, so d is on original uV scale

    es_list <- list()
    es_idx <- 1
    for (comp_idx in seq_along(whichCompareVec)) {
      comparison <- whichCompareVec[[comp_idx]]
      compData <- dataList %>% filter(blockVec %in% comparison)
      cond_name <- as.character(unique(compData$overallBlockType))

      for (sl in unique(compData$mapStimLevel)) {
        slData <- compData %>% filter(mapStimLevel == sl)
        pre_vals <- slData$PPvec[slData$blockVec == comparison[1]]
        post_vals <- slData$PPvec[slData$blockVec == comparison[2]]
        if (length(pre_vals) >= 5 & length(post_vals) >= 5) {
          d <- cohens_d(post_vals, pre_vals)
          es_list[[es_idx]] <- data.frame(
            condition = cond_name,
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

    # --- Effect size bar plot: post vs baseline by conditioning length ---
    p_es <- ggplot(es_df, aes(x = mapStimLevel, y = effect.size, fill = condition)) +
      geom_col(position = position_dodge(width = 0.7), width = 0.6) +
      geom_errorbar(aes(ymin = CI_low, ymax = CI_high),
                    position = position_dodge(width = 0.7), width = 0.2) +
      geom_hline(yintercept = c(-0.8, -0.5, -0.2, 0, 0.2, 0.5, 0.8),
                 linetype = "dashed", color = "grey50", linewidth = 0.3) +
      labs(x = "Stimulation Level",
           y = "Cohen's d (95% CI)",
           fill = "Conditioning Protocol",
           title = "Effect Size: Post vs Baseline by Conditioning Length",
           subtitle = "Computed on original EP magnitude (uV)") +
      scale_fill_brewer(palette = "Dark2") +
      theme_bw(base_size = 14)
    print(p_es)

    if (savePlot) {
      ggsave(paste0("subj_a23ed_effect_size_conditioning_length.png"),
             plot = p_es, units = "in", width = figWidth, height = figHeight, dpi = 600)
      ggsave(paste0("subj_a23ed_effect_size_conditioning_length.eps"),
             plot = p_es, units = "in", width = figWidth, height = figHeight,
             device = cairo_ps, fallback_resolution = 600)
    }

    #########
    # Permutation tests (non-parametric, no model assumptions)
    # Shuffles trial block labels within each stim level

    set.seed(42)
    nPerm <- 10000

    # --- Test 1 & 2: post vs baseline within each conditioning length ---
    # Uses median (robust to outliers) as the test statistic for both observed
    # and permuted differences. Permutation framework is valid for any statistic.
    perm_results <- data.frame(condition=character(), mapStimLevel=character(),
                               obs_diff=numeric(), perm_p=numeric(), stringsAsFactors=FALSE)

    for (comp_idx in seq_along(whichCompareVec)) {
      comparison <- whichCompareVec[[comp_idx]]
      compData <- dataList %>% filter(blockVec %in% comparison)
      cond_name <- as.character(unique(compData$overallBlockType))

      for (sl in unique(compData$mapStimLevel)) {
        slData <- compData %>% filter(mapStimLevel == sl)
        pre_vals <- slData$PPvec[slData$blockVec == comparison[1]]
        post_vals <- slData$PPvec[slData$blockVec == comparison[2]]
        obs_diff <- median(post_vals) - median(pre_vals)

        block_labels <- slData$blockVec
        all_vals <- slData$PPvec
        perm_diffs <- numeric(nPerm)
        for (p in 1:nPerm) {
          shuf <- sample(block_labels)
          perm_diffs[p] <- median(all_vals[shuf == comparison[2]]) -
                           median(all_vals[shuf == comparison[1]])
        }
        p_val <- mean(abs(perm_diffs) >= abs(obs_diff))
        perm_results <- rbind(perm_results,
          data.frame(condition=cond_name, mapStimLevel=as.character(sl),
                     obs_diff=obs_diff, perm_p=p_val))
      }
    }

    cat("\n=== Permutation tests (median): post vs baseline (per condition, per stim level) ===\n")
    print(perm_results)

    # --- Test 3: difference of differences (15-min effect - 5-min effect) ---
    # Does the post-baseline change differ between 5-min and 15-min conditioning?
    comp5 <- whichCompareVec[[1]]; comp15 <- whichCompareVec[[2]]

    perm_interaction <- data.frame(mapStimLevel=character(), obs_interaction=numeric(),
                                   perm_p=numeric(), stringsAsFactors=FALSE)

    stim_levels <- intersect(
      unique(as.character((dataList %>% filter(blockVec %in% comp5))$mapStimLevel)),
      unique(as.character((dataList %>% filter(blockVec %in% comp15))$mapStimLevel)))

    for (sl in stim_levels) {
      d5 <- dataList %>% filter(blockVec %in% comp5, mapStimLevel == sl)
      d15 <- dataList %>% filter(blockVec %in% comp15, mapStimLevel == sl)

      diff5 <- median(d5$PPvec[d5$blockVec == comp5[2]]) - median(d5$PPvec[d5$blockVec == comp5[1]])
      diff15 <- median(d15$PPvec[d15$blockVec == comp15[2]]) - median(d15$PPvec[d15$blockVec == comp15[1]])
      obs_int <- diff15 - diff5

      perm_ints <- numeric(nPerm)
      for (p in 1:nPerm) {
        shuf5 <- sample(d5$blockVec)
        shuf15 <- sample(d15$blockVec)
        pdiff5 <- median(d5$PPvec[shuf5 == comp5[2]]) - median(d5$PPvec[shuf5 == comp5[1]])
        pdiff15 <- median(d15$PPvec[shuf15 == comp15[2]]) - median(d15$PPvec[shuf15 == comp15[1]])
        perm_ints[p] <- pdiff15 - pdiff5
      }
      p_val <- mean(abs(perm_ints) >= abs(obs_int))
      perm_interaction <- rbind(perm_interaction,
        data.frame(mapStimLevel=sl, obs_interaction=obs_int, perm_p=p_val))
    }

    cat("\n=== Permutation test (median): 15-min effect minus 5-min effect ===\n")
    print(perm_interaction)

    # --- Plot: permutation null distribution for the interaction ---
    # Use the pooled (across stim levels) interaction
    all_d5 <- dataList %>% filter(blockVec %in% comp5)
    all_d15 <- dataList %>% filter(blockVec %in% comp15)
    obs_diff5 <- median(all_d5$PPvec[all_d5$blockVec==comp5[2]]) - median(all_d5$PPvec[all_d5$blockVec==comp5[1]])
    obs_diff15 <- median(all_d15$PPvec[all_d15$blockVec==comp15[2]]) - median(all_d15$PPvec[all_d15$blockVec==comp15[1]])
    obs_pooled_int <- obs_diff15 - obs_diff5

    perm_pooled <- numeric(nPerm)
    for (p in 1:nPerm) {
      shuf5 <- sample(all_d5$blockVec)
      shuf15 <- sample(all_d15$blockVec)
      pd5 <- median(all_d5$PPvec[shuf5==comp5[2]]) - median(all_d5$PPvec[shuf5==comp5[1]])
      pd15 <- median(all_d15$PPvec[shuf15==comp15[2]]) - median(all_d15$PPvec[shuf15==comp15[1]])
      perm_pooled[p] <- pd15 - pd5
    }
    pooled_p <- mean(abs(perm_pooled) >= abs(obs_pooled_int))

    null_df <- data.frame(value = perm_pooled)
    p_perm <- ggplot(null_df, aes(x = value)) +
      geom_histogram(aes(fill = abs(value) >= abs(obs_pooled_int)),
                     bins = 60, color = "grey50", show.legend = FALSE) +
      scale_fill_manual(values = c("FALSE" = "grey70", "TRUE" = "#C93312")) +
      geom_vline(xintercept = obs_pooled_int, color = "#C93312", linewidth = 1.2) +
      geom_vline(xintercept = -obs_pooled_int, color = "#C93312", linewidth = 1.2, linetype = "dashed") +
      annotate("text", x = obs_pooled_int, y = Inf, vjust = 2, hjust = -0.1,
               label = sprintf("observed = %.1f uV\np = %.4f", obs_pooled_int, pooled_p),
               color = "#C93312", size = 4) +
      labs(x = expression(paste(Delta, " (15-min effect) - (5-min effect)  [",mu,"V]")),
           y = "Count",
           title = "Permutation Test: 15-min vs 5-min Conditioning Effect",
           subtitle = sprintf("%s permutations (two-sided)", formatC(nPerm, format="d", big.mark=","))) +
      theme_bw(base_size = 14)
    print(p_perm)

    if (savePlot) {
      ggsave(paste0("subj_a23ed_perm_interaction.png"),
             plot = p_perm, units = "in", width = 7, height = 4, dpi = 600)
      ggsave(paste0("subj_a23ed_perm_interaction.eps"),
             plot = p_perm, units = "in", width = 7, height = 4,
             device = cairo_ps, fallback_resolution = 600)
    }

    # --- Per-condition permutation null distribution histograms ---
    for (comp_idx in seq_along(whichCompareVec)) {
      comparison <- whichCompareVec[[comp_idx]]
      compData <- dataList %>% filter(blockVec %in% comparison)
      cond_name <- as.character(unique(compData$overallBlockType))

      obs_diff <- median(compData$PPvec[compData$blockVec == comparison[2]]) -
                  median(compData$PPvec[compData$blockVec == comparison[1]])

      perm_diffs <- numeric(nPerm)
      for (p in 1:nPerm) {
        shuf <- sample(compData$blockVec)
        perm_diffs[p] <- median(compData$PPvec[shuf == comparison[2]]) -
                         median(compData$PPvec[shuf == comparison[1]])
      }
      p_val <- mean(abs(perm_diffs) >= abs(obs_diff))

      null_df_cond <- data.frame(value = perm_diffs)
      p_perm_cond <- ggplot(null_df_cond, aes(x = value)) +
        geom_histogram(aes(fill = abs(value) >= abs(obs_diff)),
                       bins = 60, color = "grey50", show.legend = FALSE) +
        scale_fill_manual(values = c("FALSE" = "grey70", "TRUE" = "#C93312")) +
        geom_vline(xintercept = obs_diff, color = "#C93312", linewidth = 1.2) +
        geom_vline(xintercept = -obs_diff, color = "#C93312", linewidth = 1.2, linetype = "dashed") +
        annotate("text", x = obs_diff, y = Inf, vjust = 2, hjust = -0.1,
                 label = sprintf("observed = %.1f uV\np = %.4f", obs_diff, p_val),
                 color = "#C93312", size = 4) +
        labs(x = expression(paste("Post - Baseline  [",mu,"V]")),
             y = "Count",
             title = paste0("Permutation Test: ", cond_name),
             subtitle = sprintf("%s permutations (two-sided)", formatC(nPerm, format="d", big.mark=","))) +
        theme_bw(base_size = 14)
      print(p_perm_cond)

      if (savePlot) {
        fname <- gsub("[/ ]", "_", tolower(cond_name))
        ggsave(paste0("subj_a23ed_perm_", fname, ".png"),
               plot = p_perm_cond, units = "in", width = 7, height = 4, dpi = 600)
        ggsave(paste0("subj_a23ed_perm_", fname, ".eps"),
               plot = p_perm_cond, units = "in", width = 7, height = 4,
               device = cairo_ps, fallback_resolution = 600)
      }
    }

    # --- Two-panel subplot: each comparison with its own baseline ---
    panel_data <- dataList %>%
      mutate(panel = overallBlockType,
             condition = ifelse(pre_post == "pre", "Baseline", as.character(overallBlockType)))

    p_panels <- ggplot(panel_data, aes(x = as.numeric(as.character(mapStimLevel)),
                                        y = PPvec, color = pre_post)) +
      geom_point(position = position_jitterdodge(dodge.width = 0.3), alpha = 0.6) +
      geom_smooth(method = lm, se = TRUE) +
      facet_wrap(~overallBlockType) +
      scale_color_manual(values = c("pre" = "grey50", "post" = "#D95F02"),
                         labels = c("pre" = "Baseline", "post" = "Post-conditioning")) +
      labs(x = "Stimulation Level",
           y = expression(paste("EP Magnitude (",mu,"V)")),
           color = "",
           title = "Pre vs Post Conditioning by Protocol",
           subtitle = "Each panel shows its own baseline (separate blocks)") +
      theme_bw(base_size = 14) +
      theme(legend.position = "bottom")
    print(p_panels)

    if (savePlot) {
      ggsave(paste0("subj_a23ed_two_panel_prepost.png"),
             plot = p_panels, units = "in", width = 10, height = 5, dpi = 600)
      ggsave(paste0("subj_a23ed_two_panel_prepost.eps"),
             plot = p_panels, units = "in", width = 10, height = 5,
             device = cairo_ps, fallback_resolution = 600)
    }

  }
}

setwd(oldwd)
