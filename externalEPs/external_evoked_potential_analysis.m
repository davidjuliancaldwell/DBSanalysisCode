% Script to compare conditions per subject
close all;clear all;clc
Z_ConstantsDBS_externalEPs

% subjects with good response -
%1dd75
%63ce7
%80301
%B305e

SIDS = {'80301','50ad9','b26b7','1dd75','b305e','c1c8c'};
SIDSwithout6 = {'b26b7'};

plotIt = 1;
saveIt = 1;

SIDS = {'c1c8c'};

SIDS
for sid = SIDS
    sid = sid{:};
    fprintf(['subject ' sid '\n']);
    
    if strcmp(SIDSwithout6,sid)
        maxFile = 2;
    else
        maxFile = 6;
    end
    
    for file = 1:maxFile
        fprintf(['file '  num2str(file) '\n']);
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
            case '1dd75'
                switch file
                    case 1
                        load(fullfile(OUTPUT_DIR,[sid '_stim_L_singleDBS_0_1.mat']));
                    case 2
                        load(fullfile(OUTPUT_DIR,[sid '_stim_L_singleDBS_1_2.mat']));
                    case 3
                        load(fullfile(OUTPUT_DIR,[sid '_stim_L_singleDBS_2_3.mat']));
                    case 4
                        load(fullfile(OUTPUT_DIR,[sid '_stim_L_singleDBS_3_2.mat']));
                    case 5 % another really noisy trial, cond 4 - trial 4,5 , a little strange through 8
                        load(fullfile(OUTPUT_DIR,[sid '_stim_L_singleDBS_2_1.mat']));
                    case 6 % strange artifact on one of the recordings
                        load(fullfile(OUTPUT_DIR,[sid '_stim_L_singleDBS_1_0.mat']));
                end
            case 'b305e'
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
            case 'c1c8c'
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
        end
        %%
        
        % scale by 4
        ECoG_sep = cellfun(@(x) x*4,ECoG_sep,'un',0);
        DBS_sep = cellfun(@(x) x*4,DBS_sep,'un',0);
        
        windowInt = [500 1250];
        if size(ucondition,2)>size(ucondition,1)
            ucondition = ucondition';
        end
        %plot_EPs(ucondition,DBS_sep,ECoG_sep,t,stimChans,windowInt);
        %%
        tBegin = 510;
        tEnd = 650;
        
        %%[ECoGPP,DBSPP,~,~] = extract_PP(ucondition,DBS_sep,ECoG_sep,t,stimChans,tBegin,tEnd);
        smooth = 1;
        rerefMode = '';
        channelReref = 1;
        
        if max(stimChans) > 3
            stimChans = stimChans - 4;
        end
        
        % do 1:8
        if size(ECoG_sep{1},2) == 16
            
            firstStrip = [9:16];
            secondStrip = [1:8];
            
            [ECoGPPtemp,pkLoc,trLoc] = extract_PP_peak_to_peak(ucondition,ECoG_sep,t,firstStrip,tBegin,tEnd,rerefMode,channelReref,smooth);
            ECoGPP(1:8,:) = ECoGPPtemp(1:8,:);
            pkLocVec(1:8,:) = tBegin + round(1e3.*pkLoc(1:8,:)./(ECOG_fs));
            trLocVec(1:8,:) = tBegin + round(1e3.*trLoc(1:8,:)./(ECOG_fs));
            
            [ECoGPPtemp,pkLoc,trLoc] = extract_PP_peak_to_peak(ucondition,ECoG_sep,t,secondStrip,tBegin,tEnd,rerefMode,channelReref,smooth);
            ECoGPP(9:16,:) = ECoGPPtemp(9:16,:);
            pkLocVec(9:16,:) = tBegin + round(1e3.*pkLoc(9:16,:)./(ECOG_fs));
            trLocVec(9:16,:) = tBegin + round(1e3.*trLoc(9:16,:)./(ECOG_fs));
            
        else
            % badChansVec = size(ECoG_sep{1},2);
            % badChansVec(1:8) = 0;
            badChansVec = [];
            [ECoGPPtemp,pkLoc,trLoc] = extract_PP_peak_to_peak(ucondition,ECoG_sep,t,badChansVec,tBegin,tEnd,rerefMode,channelReref,smooth);
            ECoGPP(1:8,:) = ECoGPPtemp(1:8,:);
            pkLocVec(1:8,:) = tBegin + round(1e3.*pkLoc(1:8,:)./(ECOG_fs));
            trLocVec(1:8,:) = tBegin + round(1e3.*trLoc(1:8,:)./(ECOG_fs));
        end
        
        
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % now significance
        ECoGsepbaselineAvg = {};
        DBSsepBaselineAvg = {};
        ECoGsepMean = {};
        DBSsepMean = {};
        
        % rerefMode = 'mean';
        rerefMode = '';
        for i = ucondition'
            ECoGtemp = ECoG_sep{i};
            DBStemp = DBS_sep{i};
            cellMode = {'median','bipolar','mean','bipolarPair','singleChan',''};
            
            if sum(strcmp(rerefMode,cellMode))
                if size(ECoG_sep{1},2) == 16
                    
                    firstStrip = [9:16];
                    secondStrip = [1:8];
                    
                    ECoGtemp1 = rereference_CAR_median(ECoGtemp,rerefMode,firstStrip,[1 2 3],channelReref);
                    
                    ECoGtemp2 = rereference_CAR_median(ECoGtemp,rerefMode,secondStrip,[1 2 3],channelReref);
                    
                    ECoGtemp(:,1:8,:) = ECoGtemp1(:,1:8,:);
                    ECoGtemp(:,9:16,:) = ECoGtemp2(:,9:16,:);
                    
                else
                    ECoGtemp = ECoGtemp(:,1:8,:);
                    ECoGtemp = rereference_CAR_median(ECoGtemp,rerefMode,[],[1 2 3],channelReref);
                    
                end
            end
            
            ECoGsepMean{i} = squeeze(mean(ECoGtemp,3));
            DBSsepMean{i} = squeeze(mean(DBStemp,3));
            
            ECoGsepbaselineAvg{i} = squeeze(mean(ECoGtemp(t<0,:,:),3));
            DBSsepBaselineAvg{i} = squeeze(mean(DBStemp(t<0,:,:),3));
            
        end
        
        ECoGPPzScore = ECoGPP;
        
        for i = ucondition'
            ECoGPPzScore(:,i) = normalize_pp(ECoGPP(:,i),ECoGsepbaselineAvg{i}');
        end
        %%  Plotting section raw Pk-Pk
        if plotIt
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            figure
            hold on
            elecs = 1:size(ECoGPP,1);
            for i = 1:size(ECoGPP,2)
                scatter(elecs,1e6*ECoGPP(:,i),'filled')
            end
            legend({'1','2','3','4'},'location','northwest');
            xlabel('electrodes')
            ylabel('peak to peak voltage (\muV)')
            set(gca,'fontsize',14')
            title([sid ' stimulation channels ' num2str(stimChans)])
            SaveFig(OUTPUT_DIR, [sid '_PkPk_' num2str(stimChans(1)) '_' num2str(stimChans(2))], 'png', '-r600');
            
            
            %% Plotting Z-score Pk-Pk
            
            figure
            hold on
            elecs = 1:size(ECoGPPzScore,1);
            for i = 1:size(ECoGPPzScore,2)
                scatter(elecs,ECoGPPzScore(:,i),'filled')
            end
            legend({'1','2','3','4'},'location','northwest');
            xlabel('electrodes')
            ylabel('peak to peak z-score')
            set(gca,'fontsize',14')
            title({[sid ' stimulation channels ' num2str(stimChans)], ' z score relative to baseline'})
            
            SaveFig(OUTPUT_DIR, [sid '_zScorePkPk_' num2str(stimChans(1)) '_' num2str(stimChans(2))], 'png', '-r600');
            
            %% Plotting raw time series for channel with greatest Pk-Pk in most intense stimulation condition
            [~,chanInt] = max(ECoGPP(:,4)) ;
            figure
            plot(t,1e6*ECoGsepMean{1}(:,chanInt))
            hold on
            plot(t,1e6*ECoGsepMean{2}(:,chanInt))
            plot(t,1e6*ECoGsepMean{3}(:,chanInt))
            plot(t,1e6*ECoGsepMean{4}(:,chanInt))
            
            title({[sid ' channel ' num2str(chanInt)], ' raw signal traces ', ['stimulation channels ' num2str(stimChans)]})
            legend({'1','2','3','4'});
            xlabel('time (ms)')
            ylabel('voltage (\muV)')
            % xlim([450 800])
            xlim([-100 800])
            ylim([-240 240])
            vline(pkLocVec(chanInt,4))
            vline(trLocVec(chanInt,4))
            set(gca,'fontsize',14)
            if saveIt
                SaveFig(OUTPUT_DIR, [sid '_rawTrace_stimChans_' num2str(stimChans(1)) '_' num2str(stimChans(2)) '_recordChanMax_' num2str(chanInt)], 'png', '-r600');
            end
            
            if size(ECoGtemp,2) > 8
                for chanInt = 1:16
                    figure
                    plot(t,1e6*ECoGsepMean{1}(:,chanInt))
                    hold on
                    plot(t,1e6*ECoGsepMean{2}(:,chanInt))
                    plot(t,1e6*ECoGsepMean{3}(:,chanInt))
                    plot(t,1e6*ECoGsepMean{4}(:,chanInt))
                    
                    title({[sid ' channel ' num2str(chanInt)], ' raw signal traces ', ['stimulation channels ' num2str(stimChans)]})
                    legend({'1','2','3','4'});
                    xlabel('time (ms)')
                    ylabel('voltage (\muV)')
                    %  xlim([450 800])
                    xlim([-100 800])
                    ylim([-240 240])
                    vline(pkLocVec(chanInt,4))
                    vline(trLocVec(chanInt,4))
                    set(gca,'fontsize',14)
                    if saveIt
                        SaveFig(OUTPUT_DIR, [sid '_rawTrace_stimChans_' num2str(stimChans(1)) '_' num2str(stimChans(2)) '_recordChan_' num2str(chanInt)], 'png', '-r600');
                    end
                end
            else
                for chanInt = 1:8
                    
                    figure
                    plot(t,1e6*ECoGsepMean{1}(:,chanInt))
                    hold on
                    plot(t,1e6*ECoGsepMean{2}(:,chanInt))
                    plot(t,1e6*ECoGsepMean{3}(:,chanInt))
                    plot(t,1e6*ECoGsepMean{4}(:,chanInt))
                    
                    title({[sid ' channel ' num2str(chanInt)], ' raw signal traces ', ['stimulation channels ' num2str(stimChans)]})
                    legend({'1','2','3','4'});
                    xlabel('time (ms)')
                    ylabel('voltage (\muV)')
                    % xlim([450 800])
                    xlim([-100 800])
                    ylim([-240 240])
                    vline(pkLocVec(chanInt,4))
                    vline(trLocVec(chanInt,4))
                    set(gca,'fontsize',14)
                    if saveIt
                        SaveFig(OUTPUT_DIR, [sid '_rawTrace_stimChans_' num2str(stimChans(1)) '_' num2str(stimChans(2)) '_recordChan_' num2str(chanInt)], 'png', '-r600');
                    end
                end
            end
        end
        % clear variables except ones required for loop references
        clearvars -except META_DIR OUTPUT_DIR SIDS sid plotIt file maxFile SIDSwithout6 saveIt
        close all
    end
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