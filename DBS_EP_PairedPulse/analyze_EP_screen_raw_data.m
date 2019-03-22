% script to analyze EP screen from paired pulse DBS experiments
close all;clear all;clc
Z_ConstantsDBS_PairedPulse;
%% load in subject

% here are subjects that we have acquired MEP data on
%sid = '3809e';
sid = '68574';
sid = 'e9c9b';
sid = '46c2a';
sidCell = {'46c2a','68574','e9c9b','a23ed'};
blockCell = {[1:6],[5,6,13,14],[1:4,7],[3,4,9]};
indexVec = [1,2,3,4];
indexVec = [4];
matlab_dir = 'MATLAB_Converted';
experiment = 'EP_Screen';

for index = indexVec
    
    sid = sidCell{index}
    blockVec = blockCell{index}
    for block = blockVec
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
                        
                        load(fullfile(SUB_DIR,sid,matlab_dir,experiment,'EP_Screen-2.mat'));
                        stimChans = [4 5];
                    case 3
                        load(fullfile(SUB_DIR,sid,matlab_dir,experiment,'EP_Screen-3.mat'));
                        stimChans = [4 5];
                    case 4
                        load(fullfile(SUB_DIR,sid,matlab_dir,experiment,'EP_Screen-4.mat')); % most promising one
                        stimChans = [7 8];
                    case 5
                        load(fullfile(SUB_DIR,sid,matlab_dir,experiment,'EP_Screen-5.mat'));
                        stimChans = [8 7];
                    case 6
                        load(fullfile(SUB_DIR,sid,matlab_dir,experiment,'EP_Screen-6.mat')); % most promising one
                        stimChans = [6 5];
                end
            case '68574'
                switch block
                    case 5
                        load(fullfile(SUB_DIR,sid,matlab_dir,experiment,'EP_Screen-5.mat'));
                        stimChans = [6 5];
                    case 6
                        load(fullfile(SUB_DIR,sid,matlab_dir,experiment,'EP_Screen-6.mat')); %
                        stimChans = [5 6];
                    case 13
                        load(fullfile(SUB_DIR,sid,matlab_dir,experiment,'EP_Screen-13.mat'));
                        stimChans = [6 5];
                    case 14
                        load(fullfile(SUB_DIR,sid,matlab_dir,experiment,'EP_Screen-14.mat')); %
                        stimChans = [5 6];
                end
            case 'e9c9b'
                switch block
                    case 1
                        load(fullfile(SUB_DIR,sid,matlab_dir,experiment,'EP_Screen-1.mat'));
                        stimChans = [8 7];
                    case 2
                        load(fullfile(SUB_DIR,sid,matlab_dir,experiment,'EP_Screen-2.mat')); %
                        stimChans = [7 8];
                    case 3
                        load(fullfile(SUB_DIR,sid,matlab_dir,experiment,'EP_Screen-3.mat'));
                        stimChans = [7 6];
                    case 4
                        load(fullfile(SUB_DIR,sid,matlab_dir,experiment,'EP_Screen-4.mat')); %
                        stimChans = [6 7];
                    case 7
                        load(fullfile(SUB_DIR,sid,matlab_dir,experiment,'EP_Screen-7.mat')); %
                        stimChans = [7 8];
                end
                
            case 'a23ed'
                switch block
                    case 1
                        load(fullfile(SUB_DIR,sid,matlab_dir,experiment,'EP_Screen-1.mat'));
                        stimChans = [7 8];
                    case 2
                        load(fullfile(SUB_DIR,sid,matlab_dir,experiment,'EP_Screen-2.mat')); %
                        stimChans = [8 7];
                    case 3
                        load(fullfile(SUB_DIR,sid,matlab_dir,experiment,'EP_Screen-3.mat'));
                        stimChans = [7 6];
                    case 4
                        load(fullfile(SUB_DIR,sid,matlab_dir,experiment,'EP_Screen-4.mat')); %
                        stimChans = [6 7];
                    case 5
                        load(fullfile(SUB_DIR,sid,matlab_dir,experiment,'EP_Screen-5.mat')); %
                        stimChans = [6 5];
                    case 6
                        load(fullfile(SUB_DIR,sid,matlab_dir,experiment,'EP_Screen-6.mat')); %
                        stimChans = [5 6];
                    case 7
                        load(fullfile(SUB_DIR,sid,matlab_dir,experiment,'EP_Screen-7.mat')); %
                        stimChans = [5 4];
                    case 8
                        load(fullfile(SUB_DIR,sid,matlab_dir,experiment,'EP_Screen-8.mat')); %
                        stimChans = [4 5];
                    case 9
                        load(fullfile(SUB_DIR,sid,matlab_dir,experiment,'EP_Screen-9.mat')); %
                        stimChans = [6 7];
                        
                end
        end
        
        %%
        ECoG = 4.*Wav1.data;
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
        plotIt = 0;
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
        
        preTime = 10;
        postTime = 100;
        
        preSamps = round(preTime*blckedDataFs/1e3);
        postSamps = round(postTime*blckedDataFs/1e3);
        tEpoch = [-preSamps:postSamps-1]/blckedDataFs*1e3;
        
        epochUnNormed = getEpochSignal(ECoG,sts-preSamps,sts+postSamps);
        epochNormed =  epochUnNormed-repmat(mean(epochUnNormed(tEpoch<0,:,:),1), [size(epochUnNormed, 1), 1]);
        
        epochNormed(:,stimChans,:) = 0;
        
        xlims = [-5 50];
        
        for i = uniquePulseWidthLabels
            if i(1)>1000
                numelements = size(epochNormed(:,:,labels==i(1) & pulseWidths == i(2)),3);
                fprintf(['Subject ' sid ' ' num2str(i(1)) ' ',num2str(i(2)) ' ' num2str(stimChans(1)) ' ' num2str(stimChans(2)) ' # avg ' num2str(numelements) '\n'])
                smallMultiples_DBS(epochNormed(:,:,labels==i(1) & pulseWidths == i(2)),tEpoch/1e3,'type2',stimChans,'xlims',xlims);
                sgtitle(['stim chans = ' num2str(stimChans) ' current level = ' num2str(i(1)) ' \muA pulse width = ' num2str(i(2)) ' \mus'])
            end
        end
        
    end
end