%% iterate through blocks
for block = blocks
    
    
    load(fullfile(fileFolder,sid,matlabFolder,dataFolder,['EP_Measure-' num2str(block) '.mat']));
    
    
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
    preSamps = round(50*ECoGfs/1e3);
    postSamps = round(500*ECoGfs/1e3);
  % postSamps = round(1000*ECoGfs/1e3);
    tEpoch = [-preSamps:postSamps-1]/ECoGfs*1e3;
    
    goodVec = logical(ones(size(ECoG,2),1));
    stimChans = stimChansVec(blockCount,:);
    goodVec(stimChans) = 0;
    chansList = [1:8];
    chans = chansList(goodVec);
    ECoG(:,stimChans,:) = 0;
    
    count = 1;
    
    for i = stimLevelUniq
        
        epochUnNormed = getEpochSignal(ECoG,stimLevelCell{count}-preSamps,stimLevelCell{count}+postSamps);
        epochNormed =  epochUnNormed-repmat(mean(epochUnNormed(tEpoch<0,:,:),1), [size(epochUnNormed, 1), 1]);
        epochsEP{count} = epochNormed;
        if plotCondAvg
            smallMultiples_DBS(mean(epochsEP{count},3),tEpoch/1e3,'type2',stimChans);
            title(['stim pair = '  num2str(stimChans) ' Stim Level ' num2str(i)])
            
            if savePlot
                SaveFig(OUTPUT_DIR, [sid '_EP_block_' num2str(block) '_stimChans_' num2str(stimChans(1)) '_' num2str(stimChans(2)) '_stimLev_' num2str(i)], 'png', '-r600');
            end
        end
        count = count + 1;
        
    end
        
    epochsEPblock{blockCount} = epochsEP;
    
    % do not need to extract peak to peak values in OR 
    % renable to get PP values
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % get peak to peak values
%     rerefMode = 'none';
%     smooth = 1;
%     
%     [signalPP,pkLocs,trLocs] =  extract_PP_peak_to_peak(stimLevelUniq,epochsEP,tEpoch,stimChans,tBegin,tEnd,rerefMode,[],smooth);
%     
%     signalPPblock{blockCount} = signalPP;
%     pkLocsBlock{blockCount} = pkLocs;
%     trLocsBlock{blockCount} = trLocs;
%     
%     [signalPP,pkLocs,trLocs] =  extract_PP_peak_to_peak_single_trial(stimLevelUniq,epochsEP,tEpoch,stimChans,tBegin,tEnd,rerefMode,[],smooth);
%     
%     signalPPblockST{blockCount} = signalPP;
%     pkLocsBlockST{blockCount} = pkLocs;
%     trLocsBlockST{blockCount} = trLocs;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    blockCount = blockCount + 1;
    
    clearvars badTrials badTrialLocations
    
    fprintf(['finished block '  num2str(block) '\n'])
    
end


