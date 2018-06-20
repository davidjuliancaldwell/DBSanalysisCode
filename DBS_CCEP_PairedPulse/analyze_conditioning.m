% script to analyze EP results from paired pulse DBS experiments
%close all;clear all;clc
Z_ConstantsDBS_PairedPulse;
%% load in subject
% here are subjects that we have acquired MEP data on

sid = '3809e';
block = 1;
% load in tank
switch sid
    
    case '3809e'
        switch block
            case 1
                load(fullfile(SUB_DIR,sid,'Matlab_conversions\PairedPulseConditioning\PairedPulse-1.mat'));
                stimChans = [6 5];
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
%%
figure
plot(stimBox(:,1));
figure
plot(stimBox(:,4));

%%
figure
plot(stimProgrammed(:,1));

figure
plot(stimProgrammed(:,4));