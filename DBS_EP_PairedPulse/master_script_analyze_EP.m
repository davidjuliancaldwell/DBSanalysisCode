% script to analyze EP results from paired pulse DBS experiments
close all;clear all;clc
Z_ConstantsDBS_PairedPulse;
%% load in subject
matlab_dir = 'MATLAB_Converted';
experiment = 'EP_Measurement';

avgTrialsVec = [0,1]';
avgTrialsVec = [0]';
numAvg = 3;

savePlot = 0;
saveData = 0;
screenBadChans = 0;
plotCondAvg = 0;

plotPkTr = 0;

sidVecIterate = {'46c2a','9f852','8e907','08b13'};
sidVecIterate = {'08b13'};
sidVecIterate = {'c963f'};
sidVecIterate = {'2e114'};
sidVecIterate = {'3d413'};
sidVecIterate = {'a23ed'};

for avgTrials = avgTrialsVec'
    for sid = sidVecIterate
        sid = sid{:};
        switch sid
            case '46c2a'
                blocks = [1 3];
                chanIntList = [6];
                legendText = {'baseline','post 25 ms A/B'};
                
            case 'c963f'
                blocks = [2 4 6 8];
                chanIntList = [6];
                legendText = {'baseline','post 25 ms A/B','baseline 2','post 50 ms A/B'};
                
            case '2e114'
                blocks = [1 3 5 7 10 12];
                chanIntList = [4];
                legendText = {'baseline 1','post 25 ms A/B','baseline 2','post 50 ms A/B','baseline 3','post A only'};
                
            case '3d413'
                blocks = [1 2 3 4 5 6 7 8 9 10];
                chanIntList = [4 6];
                legendText = {"asleep","asleep","awake","awake","awake","awake","asleep","asleep","asleep","asleep"};
                
            case '9f852'
                blocks = [2 3 4 5 6 7 10 11 12];
                chanIntList = [4];
                legendText = {'baseline 2 (pre conditioning)' ,'post A/B 25 ms','baseline 3 (post 25 ms)',...
                    'post A/A 25 ms','baseline 4 (post 25 ms A/A)','post A/B 200 ms','baseline 5 - post A/B 200 ms 12 minutes later','post A/B 25 ms second time','baseline 6'};
                
            case '8e907'
                blocks = [1 2 5 6 7];
                chanIntList = [4 5];
                legendText = {'baseline 1' ,'baseline 2','post A/B 200 ms 1','post A/B 200 ms 2','post A/B 200 ms 3'};
                
            case '08b13'
                blocks = [1 3 5 6 7 8];
                chanIntList = [5 6];
                legendText = {'baseline 1' ,'baseline 2','post A/B 200 ms 1',...
                    'post A/B 200 ms 2/pre A/A','post A/A 200 ms 1','post A/A 200 ms 2'};
                
            case 'e9c9b'
                blocks = [1 2 3 4 5 6 7 8 9 10];
                chanIntList = [5 6];
                legendText = {'baseline 1','baseline 2','post A/A 200 ms 1','baseline 3 (post A/A)',...
                    'post A/B 200 ms 1','baseline 4 (post A/B)','post A/B 25 ms 1','baseline 5 (post A/B)','baseline 6',...
                    'post A/A 25 ms 1'};
                
            case '41a73'
                blocks = [1 2 3 4 5 6 7 8 9 10 11 12];
                chanIntList = [5 8];
                legendText = {'baseline 1','baseline 2','post A/B 200 ms 1','baseline 3 (post A/B)','baseline 4',...
                    'post A/A 200 ms 1','baseline 5 (post A/A)','baseline 6','baseline 7','baseline 8','baseline 9','baseline 10'};
                
            case '68574'
                blocks = [1 2 3 4 5 6 7 8 9 10 11 13 14 15];
                chanIntList = [4 7];
                legendText = {'baseline 1','baseline 2','post A/A 100 ms 1','baseline 3 (post A/A)',...
                    'post A/B 100 ms 1','baseline 4 (post A/B)','post A/A 200 ms 1','baseline 5 (post A/A)','post A/B 200 ms 1',...
                    'baseline 6 (post A/B 200 ms 1)','baseline 7 - pre DBS'...
                    ,'during DBS','post DBS 1','post DBS 2'};
                
            case '01fee'
                blocks = [1 2 3 4 5 6 7 8 9 10 11];
                chanIntList = [5 8];
                legendText = {'baseline 1','baseline 2','post A/B 100 ms 1','baseline 3 (post A/B)',...
                    'post A/A 100 ms 1','baseline 4 (post A/A)','post A/B 200 ms 1','baseline 5 (post A/B, pre DBS)',...
                    'during DBS','post DBS 1','post DBS 2'};
                
            case 'a23ed'
                blocks = [1 2 3 4 5 6 7 8 9];
                chanIntList = [5 8];
                legendText = {'baseline 1','baseline 2','post A/B 200 ms 1 (5 mins)','baseline 3 (post A/B)',...
                    'baseline 4','post A/B 200 ms (15 mins)','baseline 5','post A/A 200 ms (15 mins)','baseline 6'};  
                
                         legendText = {'baseline 2','post A/B 200 ms 1 (5 mins)','post A/B 200 ms (15 mins)','post A/A 200 ms (15 mins)'};  
                
                blocks = [2 3 6 8];
                
        end
        
        fprintf([sid,'\n'])
        fprintf(['average trials ' num2str(avgTrials) '\n'])
        
        %% prepare data
        prepare_EP_blocks
        
        %% compare multiples blocks
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
                tempChannel = repmat([1:8]',[1,size(tempPP,2)]);
                tempBlock = repmat(blocks(ii),size(tempPP));
                tempStimLevel = repmat(blockLabel{ii}{jj}',[size(tempPP,1),1]);
                
                PPvec = [PPvec; tempPP(:)];
                blockVec = [blockVec; tempBlock(:)];
                stimLevelVec = [stimLevelVec; tempStimLevel(:)];
                chanVec = [chanVec; tempChannel(:)];
            end
        end
        sidVec = cellstr(repmat(sid,[size(PPvec),1]));
        %
        T = table(PPvec,blockVec,stimLevelVec,sidVec,chanVec);
        
        %% save data for statistical analysis in table form
        
        if saveData && ~avgTrials
            writetable(T,[sid '_PairedPulseData.csv'],'Delimiter',',','QuoteStrings',true)
            %    save([sid '_PairedPulseData.mat'],'signalPPblockST','chanIntList','blocks','sid','tBegin','tEnd','blockLabel','stimLevelUniq','legendText')
        elseif saveData && avgTrials
            writetable(T,[sid '_PairedPulseData_avg.csv'],'Delimiter',',','QuoteStrings',true)
            %  save([sid '_PairedPulseData_avg.mat'],'signalPPblockST','chanIntList','blocks','sid','tBegin','tEnd','blockLabel','stimLevelUniq','legendText')
        end
        
        clearvars signalPPblockST
        
    end
end