% script to analyze EP results from paired pulse DBS experiments
close all;clear all;clc
Z_ConstantsDBS_PairedPulse;
%% load in subject
% here are subjects that we have acquired MEP data on

matlab_dir = 'MATLAB_Converted';
experiment = 'EP_Measurement';

% select blocks to look at
%blocks = [6 7 9 11 13 15]; % this is for the 7/8 pair fo re6f3c
%blocks = [5 10 14 17]; % this is for the 7/8 pair fo re6f3c


%sid = '3809e';
sid = '46c2a';
%sid = 'c963f';
%sid = '2e114';
%sid = '3d413';
%sid = 'fe7df';
%sid = 'e6f3c';
%sid = '8e907';

switch sid
    case '46c2a'
        blocks = [1 3];
        chanIntList = [6];
        
        legendText = {'baseline','post 25 ms A/B'};
        
    case '8e907'
        blocks = [1 2 5 6 7];
        chanIntList = [4 5];
        legendText = {'baseline 1' ,'baseline 2','post A/B 200 ms 1 ','post A/B 200 ms 2','post A/B 200 ms 3'};
        
        
end
savePlot = 0;
saveData = 1;
plotCondAvg = 0;

fprintf([sid,'\n'])

%% prepare data
prepare_EP_blocks

%% compare multiples blocks

%legendText = {'baseline' ,'post A/B 200 ms delay ','post A/A 25 ms control','post A/B 25 ms delay'};
%legendText = {'pre conditioning (2nd pre baseline)' ,'post 25 ms A/B - first time', 'post 25 ms A/A','post 200 ms A/B','post 200 ms A/B - 12  minutes later','post 25 ms A/B - second time'}; % what is the legend text

analyze_EP_compare_multiple_blocks

%% save data for statistical analysis in table form

PPvec = [];
blockVec = [];
stimLevelVec = [];
sidVec = [];
chanVec = [];
for ii = 1:size(signalPPblockST,2)
    for jj = 1:size(signalPPblockST{ii},2)
        tempPP = signalPPblockST{ii}{jj};
        tempChannel = repmat([1:8]',1,size(tempPP,2));
        tempBlock = repmat(blocks(ii),size(tempPP));
        tempStimLevel = repmat(blockLabel{ii}{jj}',size(tempPP,1),1);
        
        PPvec = [PPvec; tempPP(:)];
        blockVec = [blockVec; tempBlock(:)];
        stimLevelVec = [stimLevelVec; tempStimLevel(:)];
        chanVec = [chanVec; tempChannel(:)];
    end
end
sidVec = cellstr(repmat(sid,size(PPvec),1));
%
T = table(PPvec,blockVec,stimLevelVec,sidVec,chanVec);

%% save data for statistical analysis in table form

if saveData
    writetable(T,[sid '_PairedPulseData.csv'],'Delimiter',',','QuoteStrings',true)
    save([sid '_PairedPulseData.mat'],'signalPPblockST','chanIntList','blocks','sid','tBegin','tEnd','blockLabel','stimLevelUniq','legendText')
end