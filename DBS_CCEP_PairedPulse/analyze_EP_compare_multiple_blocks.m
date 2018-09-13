% script to analyze EP results from paired pulse DBS experiments
close all;clear all;clc
Z_ConstantsDBS_PairedPulse;
%% load in subject
% here are subjects that we have acquired MEP data on

%sid = '3809e';
%sid = '46c2a';
%sid = 'c963f';
%sid = '2e114';
sid = '3d413';
savePlot = 1;

sid


matlab_dir = 'MATLAB_Converted';
experiment = 'EP_Measurement';

blocks = [1 3 5 7 9];
%blocks = [2 4 6 8 10];

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
            tBegin = 1.2;
            tEnd = 35;
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
            
            tBegin = 2.1; % ms
            tEnd = 25; % ms
            
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
            
            tBegin = 2.1; % ms was 2.5 for 42
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
            
        case '3d413'
            tBegin = 2.1; % ms was 2.5 for 42
            tEnd = 35; % ms % was 35 for 426
            
            switch block
                % 1st baseline - patient was under/waking up
                case 1
                    stimChans = [7 8];
                    
                    % 1st baseline - patient was under/waking up
                case 2
                    stimChans = [6 5];
                    
                    % 2nd baseline, patient awake
                case 3
                    stimChans = [7 8];
                    
                    % 2nd baseline, patient awake
                    % pain meds added
                case 4
                    stimChans = [6 5];
                    
                    % 3nd baseline, patient awake, during CT scan
                case 5
                    stimChans = [7 8];
                    
                    % 3nd baseline, patient awake
                case 6
                    stimChans = [6 5];
                    
                    %4th baseline, going under
                case 7
                    stimChans = [7 8];
                    
                    %4th baseline, going under
                case 8
                    stimChans = [6 5];
                    
                    %5th baseline, under
                case 9
                    stimChans = [7 8];
                    
                    %5th baseline, under
                case 10
                    stimChans = [6 5];
                    
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
    
    
    fac = ECoGfs/tactFs;
    
    % adjust for when no stimuli were delivered when the stim level was
    % changed
    count = 1;
    % make a copy
    stimLevelUniqAll = stimLevelUniq;
    % get rid of any that are empty
    for i = stimLevelUniqAll
        if sum(stimCommandTimes(stimLevelCommandTimes == i)) == 0
            stimLevelUniq(count) = nan;
        end
        count = count + 1;
        
    end
    
    stimLevelUniq = stimLevelUniq(~isnan(stimLevelUniq));
    
    % seems to be 29 sample delay between stim command and when the ECoG
    % data starts to move
    % and measured peak
    count = 1;
    
    for i = stimLevelUniq
        stimLevelCell{count} = round((29+stimCommandTimes(stimLevelCommandTimes == i))*fac);
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
        
        %         if savePlot
        %             SaveFig(OUTPUT_DIR, [sid '_EP_block_' num2str(block) '_stimChans_' num2str(stimChans(1)) '_' num2str(stimChans(2)) '_stimLev_' num2str(i)], 'png', '-r600');
        %         end
        %
        
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

% select channel of interest
chanIntList =  [3,4,7];
chanIntList = [6];

% set colormap
cmap = cbrewer('qual','Dark2',length(blocks));
for chanInt = chanIntList
    for condInt = 1:4
        % confidence interval
        type = 'CI';
        figure
        set(gcf,'position',[1.4277e+03 556.3333 786.6667 678.6667])
        for i = 1:length(blocks)
            plotBTLError(tEpoch, 1e6*squeeze(epochsEPblock{i}{condInt}(:,chanInt,:)),type,cmap(i,:)')
        end
        xlim([-10 50])
        ylim([-350 350])
        ylabel('Voltage (\muV)')
        xlabel('time (ms)')
        h = findobj(gca,'Type','line');
        legend([h(5),h(4),h(3),h(2),h(1)],{'under anesthesia','awake','awake','under anesthesia','under anesthesia'})
        title(['baseline variability, Channel = ' num2str(chanInt) ' , test voltage = ' num2str(stimLevelUniq(condInt)) ' \muA'])
        set(gca,'fontsize',14)
        
        
        if savePlot
            SaveFig(OUTPUT_DIR, [sid '_confInt_EP_chanRecord' num2str(chanInt) '_stimChans_' num2str(stimChans(1)) '_' num2str(stimChans(2)) '_stimLev_' num2str(stimLevelUniq(condInt))], 'png', '-r600');
        end
        %%
        % mean
        figure
        set(gcf,'position',[1.4277e+03 556.3333 786.6667 678.6667])
        
        for i = 1:length(blocks)
            plot(tEpoch, 1e6*mean(squeeze(epochsEPblock{i}{condInt}(:,chanInt,:)),2),'linewidth',2,'color',cmap(i,:)')
            hold on
        end
        xlim([-10 50])
        ylim([-350 350])
        ylabel('Voltage (\muV)')
        xlabel('time (ms)')
        %legend('pre','','post','')
        h = findobj(gca,'Type','line');
        legend([h(5),h(4),h(3),h(2),h(1)],{'under anesthesia','awake','awake','under anesthesia','under anesthesia'})
        title(['baseline variability, Channel = ' num2str(chanInt) ' , test voltage = ' num2str(stimLevelUniq(condInt)) ' \muA'])
        set(gca,'fontsize',14)
        
        
        if savePlot
            SaveFig(OUTPUT_DIR, [sid '_avg_EP_chanRecord' num2str(chanInt) '_stimChans_' num2str(stimChans(1)) '_' num2str(stimChans(2)) '_stimLev_' num2str(stimLevelUniq(condInt))], 'png', '-r600');
        end
        
        
    end
end