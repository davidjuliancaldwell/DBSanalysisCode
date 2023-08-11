% script to analyze EP results from paired pulse DBS experiments
close all;clear all;clc
Z_ConstantsDBS_PairedPulse;
%% load in subject
matlab_dir = 'MATLAB_Converted';
experiment = 'EP_Measurement';

avgTrialsVec = [0]';
%avgTrialsVec = [0 1]';
numAvg = 3;

savePlot = 1;
saveData = 0;
screenBadChans = 0;
plotCondAvg = 0;

plotPkTr = 0;

tryArtifact = 0;

sidVecIterate = {'46c2a','9f852','8e907','08b13'};
sidVecIterate = {'08b13'};
sidVecIterate = {'c963f'};
% sidVecIterate = {'2e114'};
% sidVecIterate = {'3d413'};
% sidVecIterate = {'a23ed'};
% sidVecIterate = {'fe7df','e6f3c'};
% sidVecIterate = {'e6f3c'};
% sidVecIterate = {'46c2a','9f852','8e907','08b13'};
% sidVecIterate = {'9f852'};
% 
sidVecIterate = {'46c2a','c963f','2e114','3d413','fe7df','e6f3c',...
    '9f852','8e907','08b13','e9c9b','41a73','68574',...
    '01fee','a23ed'};


%sidVecIterate = {'68574'}
%sidVecIterate = {'3d413'}

%sidVecIterate = {'46c2a'};
%sidVecIterate = {'a23ed'};

% sidVecIterate = {'2e114','3d413','fe7df','e6f3c',...
%    '9f852','8e907','08b13','e9c9b','41a73','68574',...
%    '01fee','a23ed'};


% rereference against the first channel in the array if desired
chanReref = 1;
rerefMode = 'none';


%% try artifact processing
if tryArtifact
    type = 'dictionary';
    useProcrustes = 1; % use procrustes distance for scaling artifact on each trial
    trainDuration = [0 2]; % this is how long the stimulation train was
    minDuration = 2; % minimum duration of artifact in ms
    fixedDistance = 4;
    plotIt = 0;
    useFixedEnd = 0;
    recoverExp = 0;
    
    % parameters for detecting artifact onset and offset
    pre = 0.8; % default time window to extend before the artifact pulse to ensure the artifact is appropriately detected (0.8 ms as default)
    post = 1.5; % default time window to extend before the artifact pulse to ensure the artifact is appropriately detected (1 ms as default)
    onsetThreshold = 1.5; %This value is used as absolute valued z-score threshold to determine the onset of artifacts within a train. The differentiated smoothed signal is used to determine artifact onset. This is also used in determining how many stimulation pulses are within a train, by ensuring that beginning of each artifact is within a separate artifact pulse.
    threshVoltageCut = 80; %This is used to help determine the end of each individual artifact pulse dynamically. More specifically, this is a percentile value, against which the absolute valued, z-scored smoothed raw signal is compared to find the last value which exceeds the specified percentile voltage value. Higher values of this (i.e. closer to 100) result in a shorter duration artifact, while lower values result in a longer calculated duration of artifact. This parameter therefore should likely be set higher for more transient artifacts and lower for longer artifacts.
    threshDiffCut = 80; %This is used to help determine the end of each individual artifact pulse dynamically. More specifically, this is a percentile value, against which the absolute valued, z-scored differentiated smoothed raw signal is compared to find the last value which exceeds the specified percentile voltage value. Higher values of this (i.e. closer to 100) result in a shorter duration artifact, while lower values result in a longer calculated duration of artifact. This parameter therefore should likely be set higher for more transient artifacts and lower for longer artifacts.
    
    % these are the metrics used if the dictionary method is selected. The
    % options are 'eucl', 'cosine', 'corr', for either euclidean distance,
    % cosine similarity, or correlation for clustering and template matching.
    distanceMetricDbscan = 'eucl';
    distanceMetricSigMatch = 'corr';
    
    amntPreAverage = 5; % number of samples at the beginning of each artifact pulse to use as a baseline normalization
    normalize = 'preAverage'; % method to use for normalization of each artifact pulse. 'preAverage' uses the average across the number of samples specified above
    
    % additional HDBSCAN parameters and window selection
    bracketRange = [-5:5]; %This variable sets the number of samples around the maximum voltage deflection to use for template clustering and subsequent matching. The smaller this range, the lower the dimensionality used for clustering, and the fewer points used to calculate the best matching template. This value is used to try and ensure that non-informative points are not included in the clustering and template matching. This should be set to what looks like the approximate length of the artifact's largest part.
    minPts = 2;  % Defined as k in the manuscript. This is a parameter that determines how many neighbors are used for core distance calculations for each point in the artifact window. This is a parameter that determines how many neighbors are used for core distance calculations. Increasing this parameter restricts clusters to increasingly dense areas.
    minClustSize = 3; % Defined as n in the manuscript. The minimum number of clustered artifact pulses for a cluster to be labelled as a true cluster. Increasing this number can reduce the number of clusters, and merges some clusters together that would have otherwise been labelled as individual clusters.
    outlierThresh = 0.99; % Outlier parameter for labeling artifact pulses as noise in the HDBSCAN clustering. Any artifact pulse with an outlier score greater than this will be labelled as noise. Increasing this value results in fewer points being labelled as noise
    xlimsWave = [-50,400];
end

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
                legendText = {'baseline','post 25 ms A/B','baseline 2',...
                    'post 50 ms A/B'};
                
            case '2e114'
                blocks = [1 3 5 7 10 12];
                chanIntList = [4];
                legendText = {'baseline 1','post 25 ms A/B','baseline 2',...
                    'post 50 ms A/B','baseline 3','post A only'};
                
            case '3d413'
                blocks = [1 2 3 4 5 6 7 8 9 10];
                chanIntList = [4 6];
                legendText = {"asleep","asleep","awake","awake",...
                    "awake","awake","asleep","asleep","asleep","asleep"};
                
            case 'fe7df'
                blocks = [1 2 3 4 5 6 7 8 9];
                chanIntList = [6];
                legendText = {"baseline 1","baseline 2","baseline 3","post 25 ms A/B 1","baseline 4 post 25 ms A/B 2","post 200 ms A/B 1","baseline 5 post 200 ms A/B 2","baseline 6 post 200 ms A/B 3","baseline 7 post 200 ms A/B 4"};
                
            case 'e6f3c'
                blocks = [5 6 7 8 9 10 11 12 13 14 15 16 17];
                chanIntList = [7 6];
                legendText = {"baseline 1 6/5 pair","baseline 1 8/7 pair",...
                    "baseline 2 6/5 pair","baseline 2 8/7 pair","post A/B 200 ms 1 8/7 pair",...
                    "post A/B 200 ms 1 6/5 pair","baseline 3 8/7 pair","baseline 3 6/5 pair",...
                    "post A/A 25 ms 8/7 pair","post A/A 25 ms 6/5 pair","post A/B 25 ms 8/7 pair",...
                    "post A/B 25 ms 6/5 pair - noisy","post A/B 25 ms 6/ pair - noisy"};
                
            case '9f852'
                blocks = [2 3 4 5 6 7 10 11 12];
                chanIntList = [4];
                legendText = {'baseline 2 (pre conditioning)' ,'post A/B 25 ms',...
                    'baseline 3 (post 25 ms)','post A/A 25 ms',...
                    'baseline 4 (post 25 ms A/A)','post A/B 200 ms',...
                    'baseline 5 - post A/B 200 ms 12 minutes later',...
                    'post A/B 25 ms second time','baseline 6'};
                
            case '8e907'
                blocks = [1 2 5 6 7];
                chanIntList = [4 5];
                legendText = {'baseline 1' ,'baseline 2','post A/B 200 ms 1',...
                    'post A/B 200 ms 2','post A/B 200 ms 3'};
                
            case '08b13'
                blocks = [1 3 5 6 7 8];
                chanIntList = [5 6];
                legendText = {'baseline 1' ,'baseline 2','post A/B 200 ms 1',...
                    'post A/B 200 ms 2/pre A/A','post A/A 200 ms 1',...
                    'post A/A 200 ms 2'};
                
            case 'e9c9b'
                blocks = [1 2 3 4 5 6 7 8 9 10];
                chanIntList = [5 6];
                legendText = {'baseline 1','baseline 2','post A/A 200 ms 1',...
                    'baseline 3 (post A/A)','post A/B 200 ms 1',...
                    'baseline 4 (post A/B)','post A/B 25 ms 1','baseline 5 (post A/B)','baseline 6',...
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
                
                %                          legendText = {'baseline 2','post A/B 200 ms 1 (5 mins)','post A/B 200 ms (15 mins)','post A/A 200 ms (15 mins)'};
                %
                %                 blocks = [2 3 6 8];
                
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
        PPfromAvgVec = [];
        rmsVec = [];
        for ii = 1:size(signalPPblockST,2)
            for jj = 1:size(signalPPblockST{ii},2)
                tempPP = signalPPblockST{ii}{jj};
                tempChannel = repmat([1:8]',[1,size(tempPP,2)]);
                tempBlock = repmat(blocks(ii),size(tempPP));
                tempStimLevel = repmat(blockLabel{ii}{jj}',[size(tempPP,1),1]);
                tempPPfromAvgVec = signalPPblockSTfromAvg{ii}{jj};
                tempRMS = signalRMSblock{ii}{jj};
                
                PPvec = [PPvec; tempPP(:)];
                blockVec = [blockVec; tempBlock(:)];
                stimLevelVec = [stimLevelVec; tempStimLevel(:)];
                chanVec = [chanVec; tempChannel(:)];
                PPfromAvgVec = [PPfromAvgVec; tempPPfromAvgVec(:)];
                rmsVec = [rmsVec; tempRMS(:)];
            end
        end
        sidVec = cellstr(repmat(sid,[size(PPvec),1]));
        %
        T = table(PPvec,blockVec,stimLevelVec,sidVec,chanVec,PPfromAvgVec,rmsVec);
        
        %% save data for statistical analysis in table form
        
        if saveData && ~avgTrials
            writetable(T,[sid '_PairedPulseData_new_rms_pk_pk.csv'],'Delimiter',',','QuoteStrings',true)
            %    save([sid '_PairedPulseData.mat'],'signalPPblockST','chanIntList','blocks','sid','tBegin','tEnd','blockLabel','stimLevelUniq','legendText')
        elseif saveData && avgTrials
            writetable(T,[sid '_PairedPulseData_new_rms_pk_pk_avg_3.csv'],'Delimiter',',','QuoteStrings',true)
            %  save([sid '_PairedPulseData_avg.mat'],'signalPPblockST','chanIntList','blocks','sid','tBegin','tEnd','blockLabel','stimLevelUniq','legendText')
        end
        
        clearvars signalPPblockST signalPPblockSTfromAvg
        
        if savePlot
            close all
        end
        
    end
end