% script to analyze EP screen from paired pulse DBS experiments
close all;clear all;clc
Z_ConstantsDBS_PairedPulse;
%% load in subject

% here are subjects that we have acquired MEP data on

%sid = '3809e';
sid = '46c2a';
block = 1;

sid

matlab_dir = 'MATLAB_Converted';
experiment = 'EP_Screen';

for block = 1:6
    
    % load in tank
    switch sid
        case '3972f'
            load(fullfile(SUB_DIR,sid,'Converted\3972f_CCEP_converted\CCEP-1.mat'));
            stimChans = [7 8];
        case '3809e'
            switch block
                case 1 
                    load(fullfile(SUB_DIR,sid,matlab_dir,experiment,'EP_Screen-1.mat'));
                    stimChans = [8 7];
                case 2        
                    load(fullfile(SUB_DIR,sid,matlab_dir,experiment,'EP_Screen-2.mat'));
                    stimChans = [8 7];
                case 3
                    load(fullfile(SUB_DIR,sid,matlab_dir,experiment,'EP_Screen-3.mat'));
                    stimChans = [7 6];
                case 4
                    load(fullfile(SUB_DIR,sid,matlab_dir,experiment,'EP_Screen-4.mat')); % most promising one
                    stimChans = [6 5];
            end
        case '46c2a'
            switch block
                case 1
                    load(fullfile(SUB_DIR,sid,matlab_dir,experiment,'EP_Screen-1.mat'));
                    stimChans = [6 7];
                case 2
                    
                    load(fullfile(SUB_DIR,sid,matlab_dir,experiment,'EP_Screen-3.mat'));
                    stimChans = [4 5];
                case 3
                    load(fullfile(SUB_DIR,sid,matlab_dir,experiment,'EP_Screen-3.mat'));
                    stimChans = [4 5];
                case 4
                    load(fullfile(SUB_DIR,sid,matlab_dir,experiment,'EP_Screen-4.mat')); % most promising one
                    stimChans = [7 8];
                case 5
                    load(fullfile(SUB_DIR,sid,matlab_dir,experiment,'EP_Screen-3.mat'));
                    stimChans = [8 7];
                case 6
                    load(fullfile(SUB_DIR,sid,matlab_dir,experiment,'EP_Screen-4.mat')); % most promising one
                    stimChans = [6 5];
            end
    end
    
    %%
    ECoG = Wav1.data;
    ECoGfs = Wav1.info.SamplingRateHz;
    
    stimBox = Stim.data;
    stimFs = Stim.info.SamplingRateHz;
    stimProgrammed = Sing.data;
    
    
    blckedData = Blck.data;
    blckedDataFs = Blck.info.SamplingRateHz;
    preTime = 10;
    postTime = 30;
    preSamps = round(preTime*blckedDataFs/1e3);
    postSamps = round(postTime*blckedDataFs/1e3);
    tBlck = [-preSamps:postSamps-1]/blckedDataFs*1e3;
    blckedData = blckedData(:,1:8,:);
    
    %%
    plotIt = 1;
    savePlot = 0;
    [stim1Epoched,t,fs,labels,pulseWidths,uniqueLabels,uniquePulseWidths,uniquePulseWidthLabels] = voltage_monitor_different_width(Stim,Sing,plotIt,savePlot,'','','',1);
    
    %% find out which each of the programmed stimuli actually were set to be delivered
    [sts,bursts] = get_epoch_indices(Sing.data,ECoGfs,stimFs);

    %%
    goodVec = logical(ones(size(blckedData,2),1));
    goodVec(stimChans) = 0;
    chansList = [1:8];
    chans = chansList(goodVec);
    blckedData(:,stimChans,:) = 0;
    xlims = [-5 30];
    
    for i = uniquePulseWidthLabels
        smallMultiples_DBS(blckedData(:,:,labels==i(1) & pulseWidths == i(2)),tBlck/1e3,'type2',stimChans,'xlims',xlims);
        subtitle(['stim chans = ' num2str(stimChans) ' current level = ' num2str(i(1)) ' \muA pulse width = ' num2str(i(2)) ' \mus'])
    end
    
end