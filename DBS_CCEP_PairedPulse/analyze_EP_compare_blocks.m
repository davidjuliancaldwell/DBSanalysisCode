% script to analyze EP results from paired pulse DBS experiments
close all;clear all;clc
Z_ConstantsDBS_PairedPulse;
%% load in subject
% here are subjects that we have acquired MEP data on

%sid = '3809e';
%sid = '46c2a';
%sid = 'c963f';
sid = '2e114';
savePlot = 0;

sid


matlab_dir = 'MATLAB_Converted';
experiment = 'EP_Measurement';

blocks = [5 10];
blockCount = 1;
for block = blocks
    
    % load in tank
    switch sid
        
        case '3809e'
            switch block
                % pre
                case 1
                    stimChans = [6 5];
                case 2
                    stimChans = [8 7];
                    % post
                case 3
                    stimChans = [6 5];
                case 4
                    stimChans = [8 7];
            end
        case '46c2a'
            switch block
                case 1
                    stimChans = [7 8];
                case 2
                    stimChans = [6 5];
                    
                    % post
                case 3
                    stimChans = [7 8];
                case 4
                    stimChans = [6 5];
            end
        case 'c963f'
            
                tBegin = 3.3; % ms was 2.5 for 42
    tEnd = 25; % ms % was 35 for 426
    
            switch block
                % pre first time
                case 1
                    stimChans = [6 5];
                case 2
                    stimChans = [7 8];
                    badTrials = 1;
                    badTrialLocations = [1:16];
                    % post first time
                case 3
                    stimChans = [6 5];
                case 4
                    stimChans = [7 8];
                    % pre second time
                    
                case 5
                    stimChans = [6 5];
                case 6
                    stimChans = [7 8];
                    % post second time
                case 7
                    stimChans = [6 5];
                case 8
                    stimChans = [7 8];
                    
            end
        case '2e114'
            
                tBegin = 3.3; % ms was 2.5 for 42
    tEnd = 35; % ms % was 35 for 426
            
            switch block
                % 1st baseline
                case 1
                    stimChans = [5 6];
                    
                    
                    % 1st baseline
                case 2
                    stimChans = [4 3];
                    
                    % 1st post stim
                case 3
                    stimChans = [5 6];
                    
                    % 1st post stim
                case 4
                    stimChans = [4 3];
                    
                    % 2nd baseline
                case 5
                    stimChans = [5 6];
                    
                    % 2nd baseline
                case 6
                    stimChans = [4 3];
                    
                    % 2nd post stim
                case 7
                    stimChans = [5 6];
                    
                    % 2nd post stim
                case 8
                    % they were stimulating the entire time
                    stimChans = [4 3];
                    
                    % 2nd post stim
                case 9
                    % this was a second "post measurement"
                    stimChans = [4 3];
                    
                    % third baseline
                case 10
                    % left side DBS connected
                    stimChans = [5 6];
                    
                    % third baseline
                case 11
                    stimChans = [ 4 3];
                    
                    % 3rd post stim
                case 12
                    stimChans = [5 6];
                    
                    % 3rd post stim
                case 13
                    stimChans = [4 3];
                    
                    
            end
            
            load(fullfile(SUB_DIR,sid,matlab_dir,experiment,['EP_Measure-' num2str(block) '.mat']));
            
    end
    
    %%
    ECoG = 4*ECO1.data(:,1:8);
    ECoGfs = ECO1.info.SamplingRateHz;
    
    stimBox = Stim.data;
    stimFs = Stim.info.SamplingRateHz;
    stimProgrammed = Sing.data;
    
    tactFs = Tact.info.SamplingRateHz;
    stimCommand = Tact.data(:,1);
    stimLevel = Tact.data(:,2);
    
    [stimCommandTimes,~] = find(stimCommand == 1);
    stimLevelCommandTimes = stimLevel(stimCommandTimes);
    stimLevelUniq = unique(stimLevel(stimLevel>0))';
    stimLevelCell = {};
    
    if exist('badTrials','var')
        stimLevelCommandTimes(badTrialLocations) = [];
        stimCommandTimes(badTrialLocations) = [];
    end
    
    count = 1;
    
    fac = ECoGfs/tactFs;
    for i = stimLevelUniq
        stimLevelCell{count} = round(stimCommandTimes(stimLevelCommandTimes == i)*fac);
        count = count + 1;
    end
    
    %%
    epochsEP = {};
    count = 1;
    preSamps = round(50*ECoGfs/1e3);
    postSamps = round(500*ECoGfs/1e3);
    postSamps = round(1000*ECoGfs/1e3);
    tEpoch = [-preSamps:postSamps-1]/ECoGfs*1e3;
    
    goodVec = logical(ones(size(ECoG,2),1));
    goodVec(stimChans) = 0;
    chansList = [1:8];
    chans = chansList(goodVec);
    ECoG(:,stimChans,:) = 0;
    
    
    for i = stimLevelUniq
        
        epochUnNormed = getEpochSignal(ECoG,stimLevelCell{count}-preSamps,stimLevelCell{count}+postSamps);
        epochNormed =  epochUnNormed-repmat(mean(epochUnNormed(tEpoch<0,:,:),1), [size(epochUnNormed, 1), 1]);
        epochsEP{count} = epochNormed;
        smallMultiples_DBS(epochsEP{count},tEpoch/1e3,'type2',stimChans);
        count = count + 1;
        title(['stim pair = '  num2str(stimChans) ' Stim Level ' num2str(i)])
        
        if savePlot
            SaveFig(OUTPUT_DIR, [sid '_EP_block_' num2str(block) '_stimChans_' num2str(stimChans(1)) '_' num2str(stimChans(2)) '_stimLev_' num2str(i)], 'png', '-r600');
        end
        
        
    end
    
    
    epochsEPblock{blockCount} = epochsEP;
    
    
    % get peak to peak values
    rerefMode = 'none';
    smooth = 1;
    
    [signalPP,pkLocs,trLocs] =  extract_PP_peak_to_peak(stimLevelUniq,epochsEP,tEpoch,stimChans,tBegin,tEnd,rerefMode,[],smooth);
    
    signalPPblock{blockCount} = signalPP;
    pkLocsBlock{blockCount} = pkLocs;
    trLocsBlock{blockCount} = trLocs;
    
    [signalPP,pkLocs,trLocs] =  extract_PP_peak_to_peak_single_trial(stimLevelUniq,epochsEP,tEpoch,stimChans,tBegin,tEnd,rerefMode,[],smooth);
    
    signalPPblockST{blockCount} = signalPP;
    pkLocsBlockST{blockCount} = pkLocs;
    trLocsBlockST{blockCount} = trLocs;
    
    blockCount = blockCount + 1;
    
    clearvars badTrials badTrialLocations
    
    
end

%% now compare

for condInt = 1:4
    %condInt = 2;
    chanInt = 4;
    type = 'CI';
    figure
    plotBTLError(tEpoch, 1e6*squeeze(epochsEPblock{1}{condInt}(:,chanInt,:)),type,'r')
    plotBTLError(tEpoch, 1e6*squeeze(epochsEPblock{2}{condInt}(:,chanInt,:)),type,'b')
    xlim([-10 50])
    ylim([-240 240])
    ylabel('Voltage (\muV)')
    xlabel('time (ms)')
    %legend('pre','','post','')
    h = findobj(gca,'Type','line');
    
    legend([h(2),h(1)],{'second baseline','third baseline'})
    title(['second baseline vs third baseline, Channel = ' num2str(chanInt) ' , test voltage = ' num2str(stimLevelUniq(condInt)) ' \muA'])
    
    percentChange = 100*(signalPPblock{2}(chanInt,condInt) - signalPPblock{1}(chanInt,condInt))/signalPPblock{1}(chanInt,condInt);
    text(10,120,{['percent change in peak to peak'], ['amplitude = ' num2str(percentChange) ' %']},'fontsize',14)
    
    [~,p] = ttest2(signalPPblockST{1}{condInt}(chanInt,:),signalPPblockST{2}{condInt}(chanInt,:));
    text(10,180,['p value = ' num2str(p)],'fontsize',14)
    
    set(gca,'fontsize',14)
%     
%     figure
% plot(tEpoch, 1e6*mean(squeeze(epochsEPblock{1}{condInt}(:,chanInt,:)),2),'r','linewidth',2)
% hold on
% plot(tEpoch, 1e6*mean(squeeze(epochsEPblock{2}{condInt}(:,chanInt,:)),2),'b','linewidth',2)
% xlim([-10 50])
% ylim([-250 250])
% ylabel('Voltage (\muV)')
% xlabel('time (ms)')
%     title(['50 ms vs 25 ms paired pulse delay, Channel = ' num2str(chanInt) ' , test voltage = ' num2str(stimLevelUniq(condInt)) ' \muA'])
% 
%     set(gca,'fontsize',14)
%     legend('baseline','test')
end
%%

