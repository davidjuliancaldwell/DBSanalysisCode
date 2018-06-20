% Script to compare conditions per subject
%close all;clear all;clc
Z_ConstantsDBS_externalEPs

% subjects with good response -
%1dd75
%63ce7
%80301
%B305e

sid = '80301';
%file = 3;
plotIt = 1;
for file = 1:6
    
    switch sid
        case '80301'
            switch file
                case 1
                    load(fullfile(OUTPUT_DIR,[sid '_stim_R_bothDBS_5_4.mat']));
                case 2
                    load(fullfile(OUTPUT_DIR,[sid '_stim_R_bothDBS_6_5.mat']));
                case 3
                    load(fullfile(OUTPUT_DIR,[sid '_stim_R_bothDBS_7_6.mat']));
                case 4
                    load(fullfile(OUTPUT_DIR,[sid '_stim_R_bothDBS_4_5.mat']));
                case 5
                    load(fullfile(OUTPUT_DIR,[sid '_stim_R_bothDBS_5_6.mat']));
                case 6
                    load(fullfile(OUTPUT_DIR,[sid '_stim_R_bothDBS_6_7.mat']));
            end
        case '50ad9'
            switch file
                case 1
                    load(fullfile(OUTPUT_DIR,[sid '_stim_L_singleDBS_0_1.mat']));
                case 2
                    load(fullfile(OUTPUT_DIR,[sid '_stim_L_singleDBS_1_2.mat']));
                case 3
                    load(fullfile(OUTPUT_DIR,[sid '_stim_L_singleDBS_2_3.mat']));
                case 4
                    load(fullfile(OUTPUT_DIR,[sid '_stim_L_singleDBS_3_2.mat']));
                case 5
                    load(fullfile(OUTPUT_DIR,[sid '_stim_L_singleDBS_2_1.mat']));
                case 6
                    load(fullfile(OUTPUT_DIR,[sid '_stim_L_singleDBS_1_0.mat']));
            end
            %     case 'bb908'
            %         switch file
            %             case 1
            %                 load(fullfile(OUTPUT_DIR,[sid '_stim_l_bothDBS_2_3.mat']));
            %             case 2
            %         end
            
        case 'b26b7'
            switch file
                case 1
                    load(fullfile(OUTPUT_DIR,[sid '_stim_l_bothDBS_2_3.mat']));
                case 2
                    load(fullfile(OUTPUT_DIR,[sid '_stim_l_bothDBS_3_2.mat']));
            end
    end
    %%
    
    windowInt = [500 1250];
    %plot_EPs(ucondition,DBS_sep,ECoG_sep,t,stimChans,windowInt);
    %%
    tBegin = 510;
    tEnd = 575;
    
    %%[ECoGPP,DBSPP,~,~] = extract_PP(ucondition,DBS_sep,ECoG_sep,t,stimChans,tBegin,tEnd);
    smooth = 1;
    rerefMode = 'mean';
    channelReref = 1;
    
    % do 1:8
    if size(ECoG_sep{1},2) == 16
        firstStrip = logical(ones(16,1));
        secondStrip = firstStrip;
        firstStrip(1:8) = 0;
        secondStrip(9:16) = 0;
        
        [ECoGPPtemp,pkLoc,trLoc] = extract_PP_peak_to_peak(ucondition,ECoG_sep,t,firstStrip,tBegin,tEnd,rerefMode,channelReref,smooth);
        ECoGPP(1:8,:) = ECoGPPtemp(1:8,:);
        pkLocVec(1:8) = pkLoc(1:8);
        trLocVec(1:8) = trLoc(1:8);
        
        [ECoGPPtemp,pkLoc,trLoc] = extract_PP_peak_to_peak(ucondition,ECoG_sep,t,secondStrip,tBegin,tEnd,rerefMode,channelReref,smooth);
        ECoGPP(9:16,:) = ECoGPPtemp(9:16,:);
        pkLocVec(9:16) = pkLoc(9:16);
        trLocVec(9:16) = trLoc(9:16);
        
    else
        [ECoGPPtemp,pkLoc,trLoc] = extract_PP_peak_to_peak(ucondition,ECoG_sep,t,secondStrip,tBegin,tEnd,rerefMode,channelReref,smooth);
        ECoGPP(1:8,:) = ECoGPPtemp(1:8,:);
        pkLocVec(1:8) = pkLoc(1:8);
        trLocVec(1:8) = trLoc(1:8);
    end
    
    
    figure
    elecs = 1:size(ECoGPP,1);
    plot(elecs,ECoGPP,'o')
    legend({'1','2','3','4'});
    xlabel('electrodes')
    ylabel('peak to peak')
    set(gca,'fontsize',14')
    title([sid ' file ' num2str(file)])
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % now significance
    ECoGsepbaselineAvg = {};
    DBSsepBaselineAvg = {};
    ECoGsepMean = {};
    DBSsepMean = {};
    
    mode = 'mean';
    
    for i = ucondition'
        ECoGtemp = ECoG_sep{i};
        DBStemp = DBS_sep{i};
        
        cellMode = {'median','bipolar','mean','bipolarPair','singleChan'};
        
        if sum(strcmp(rerefMode,cellMode))
            tempSignal = rereference_CAR_median(tempSignal,rerefMode,badChans,[1 2],channelReref);
        end
        
        ECoGsepMean{i} = squeeze(mean(ECoGtemp,3));
        DBSsepMean{i} = squeeze(mean(DBStemp,3));
        
        ECoGsepbaselineAvg{i} = squeeze(mean(ECoGtemp(t<0,:,:),3));
        DBSsepBaselineAvg{i} = squeeze(mean(DBStemp(t<0,:,:),3));
          
    end
    
    ECoGPPzScore = ECoGPP;
    
    for i = ucondition'
        ECoGPPzScore(:,i) = normalize_plv(ECoGPP(:,i),ECoGsepbaselineAvg{i}');
    end
    
    figure
    elecs = 1:size(ECoGPPzScore,1);
    plot(elecs,ECoGPPzScore,'o')
    legend({'1','2','3','4'});
    xlabel('electrodes')
    ylabel('peak to peak')
    set(gca,'fontsize',14')
    title({[sid ' file ' num2str(file)], ' z score relative to baseline'})
    
end

% %% look at fft of DBS electrodes as confirmation
% [f,p1] = fft_compute(dbs_fs,DBS_sep{1}(t>500,:,:));
%
% p1Avg = mean(p1,3);
%
% figure
% for chan = 1:8
%     subplot(8,1,chan)
%     plot(f,p1Avg(:,chan))
%     xlim([0 200])
%     ylim([0 5e-6])
% end