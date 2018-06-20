% script to analyze EP screen from paired pulse DBS experiments
close all;clear all;clc
Z_ConstantsDBS_PairedPulse;
%% load in subject

% here are subjects that we have acquired MEP data on

sid = '3809e';
block = 2;

% load in tank
switch sid
    case '3972f'
        load(fullfile(SUB_DIR,sid,'Converted\3972f_CCEP_converted\CCEP-1.mat'));
        stimChans = [7 8];
    case '3809e'
        switch block
            case 1
                
                load(fullfile(SUB_DIR,sid,'Matlab_conversions\EP_Screen\EP_Screen-3.mat'));
                stimChans = [8 7];
            case 2
                
                load(fullfile(SUB_DIR,sid,'Matlab_conversions\EP_Screen\EP_Screen-3.mat'));
                stimChans = [8 7];
            case 3
                load(fullfile(SUB_DIR,sid,'Matlab_conversions\EP_Screen\EP_Screen-3.mat'));
                stimChans = [7 6];
            case 4
                load(fullfile(SUB_DIR,sid,'Matlab_conversions\EP_Screen\EP_Screen-4.mat')); % most promising one
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

goodVec = logical(ones(size(blckedData,2),1));
goodVec(stimChans) = 0;
chansList = [1:8];
chans = chansList(goodVec);
blckedData(:,stimChans,:) = 0;

smallMultiples(blckedData,tBlck/1e3,'type2',stimChans);