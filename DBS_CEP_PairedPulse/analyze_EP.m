% script to analyze EP results from paired pulse DBS experiments
close all;clear all;clc
Z_ConstantsDBS_PairedPulse;
%% load in subject
% here are subjects that we have acquired MEP data on

sid = '3809e';
sid = '46c2a';
savePlot = 0;

sid

block = 1;

matlab_dir = 'MATLAB_Converted';
experiment = 'EP_Measurement';

for block = 1:4
    % load in tank
    switch sid
        
             case '3809e'
            switch block
                % pre
                case 1
                    load(fullfile(SUB_DIR,sid,matlab_dir,experiment,'EP_Measure-1.mat'));
                    stimChans = [6 5];
                case 2
                    load(fullfile(SUB_DIR,sid,matlab_dir,experiment,'EP_Measure-2.mat'));
                    stimChans = [8 7];
                    % post
                case 3
                    load(fullfile(SUB_DIR,sid,matlab_dir,experiment,'EP_Measure-3.mat'));
                    stimChans = [6 5];
                case 4
                    load(fullfile(SUB_DIR,sid,matlab_dir,experiment,'EP_Measure-4.mat'));
                    stimChans = [8 7];
            end
        case '46c2a'
            switch block
                case 1
                    load(fullfile(SUB_DIR,sid,matlab_dir,experiment,'EP_Measure-1.mat'));
                    stimChans = [7 8];
                case 2
                    load(fullfile(SUB_DIR,sid,matlab_dir,experiment,'EP_Measure-2.mat'));
                    stimChans = [6 5];
                    
                    % post
                case 3
                    load(fullfile(SUB_DIR,sid,matlab_dir,experiment,'EP_Measure-3.mat'));
                    stimChans = [7 8];
                case 4
                    load(fullfile(SUB_DIR,sid,matlab_dir,experiment,'EP_Measure-4.mat'));
                    stimChans = [6 5];
            end
        case 'c963f'
            switch block
                % pre first time
                case 1
                    load(fullfile(SUB_DIR,sid,matlab_dir,experiment,'EP_Measure-1.mat'));
                    stimChans = [6 5];
                case 2
                    load(fullfile(SUB_DIR,sid,matlab_dir,experiment,'EP_Measure-2.mat'));
                    stimChans = [7 8];
                    badTrials = 1;
                    badTrialLocations = [1:16];
                    % post first time
                case 3
                    load(fullfile(SUB_DIR,sid,matlab_dir,experiment,'EP_Measure-3.mat'));
                    stimChans = [6 5];
                case 4
                    load(fullfile(SUB_DIR,sid,matlab_dir,experiment,'EP_Measure-4.mat'));
                    stimChans = [7 8];
                    % pre second time
                    
                case 5
                    load(fullfile(SUB_DIR,sid,matlab_dir,experiment,'EP_Measure-5.mat'));
                    stimChans = [6 5];
                case 6
                    load(fullfile(SUB_DIR,sid,matlab_dir,experiment,'EP_Measure-6.mat'));
                    stimChans = [7 8];
                    % post second time
                case 7
                    load(fullfile(SUB_DIR,sid,matlab_dir,experiment,'EP_Measure-7.mat'));
                    stimChans = [6 5];
                case 8
                    load(fullfile(SUB_DIR,sid,matlab_dir,experiment,'EP_Measure-8.mat'));
                    stimChans = [7 8];
                    
            end
    end
    
    %%
    ECoG = ECO1.data(:,1:8);
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
    
    count = 1;
    
    fac = ECoGfs/tactFs;
    
        
    % seems to be 29 sample delay between stim command and when the ECoG
    % data starts to move
    % and measured peak
    for i = stimLevelUniq
        stimLevelCell{count} = round((29+stimCommandTimes(stimLevelCommandTimes == i))*fac);
        count = count + 1;
    end
    
    %%
    epochsEP = {};
    epochsEPpeakToPeak = [];
    count = 1;
    preSamps = round(50*ECoGfs/1e3);
    postSamps = round(500*ECoGfs/1e3);
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
    
    % get peak to peak values
    tBegin = 2.5; % ms
    tEnd = 50; % ms
    rerefMode = 'none';
    smooth = 1;
    
    [signalPP,pkLocs,trLocs] =  extract_PP_peak_to_peak(stimLevelUniq,epochsEP,tEpoch,stimChans,tBegin,tEnd,rerefMode,[],smooth);
    
    
    
    
end