% script to analyze CCEP screen from paired pulse DBS experiments 
close all;clear all;clc
Z_ConstantsDBS_PairedPulse;
%% load in subject

% here are subjects that we have acquired MEP data on

sid = '3972f';
% load in tank
switch sid
    case '3972f'
        load(fullfile(SUB_DIR,sid,'Converted\3972f_CCEP_converted\CCEP-1.mat'));
        
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
        tBlck = [-preSamps:postSamps-1]*blckedDataFs/1e3;
        stimChans = [7 8];
        blckedData = blckedData(:,1:8,:);
        
end

%%

goodVec = logical(ones(size(blckedData,2),1));
goodVec(stimChans) = 0;
chansList = [1:8]
chans = chansList(goodVec);
blckedData(:,stimChans,:) = 0;

smallMultiples(blckedData,tBlck,'type1',chans,'type2',stimChans);