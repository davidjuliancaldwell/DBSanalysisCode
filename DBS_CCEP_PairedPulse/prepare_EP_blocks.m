%% iterate through blocks
for block = blocks
    
    % load in tank
    switch sid
        
        case '3809e'
            switch block
                % pre
                case 1
                    stimChans = [6 5];
                case 2
                    stimChans = [8 7];
                    % post
                case 3
                    stimChans = [6 5];
                case 4
                    stimChans = [8 7];
            end
        case '46c2a'
            tBegin = 1.2;
            tEnd = 35;
            switch block
                case 1
                    stimChans = [7 8];
                case 2
                    stimChans = [6 5];
                    
                    % post
                case 3
                    stimChans = [7 8];
                case 4
                    stimChans = [6 5];
            end
        case 'c963f'
            
            tBegin = 2.1; % ms
            tEnd = 25; % ms
            
            switch block
                % pre first time
                case 1
                    stimChans = [6 5];
                case 2
                    stimChans = [7 8];
                    badTrials = 1;
                    badTrialLocations = [1:16];
                    % post first time
                case 3
                    stimChans = [6 5];
                case 4
                    stimChans = [7 8];
                    % pre second time
                    
                case 5
                    stimChans = [6 5];
                case 6
                    stimChans = [7 8];
                    % post second time
                case 7
                    stimChans = [6 5];
                case 8
                    stimChans = [7 8];
                    
            end
        case '2e114'
            
            tBegin = 1.8; % ms was 2.5 for 42
            tEnd = 35; % ms % was 35 for 426
            
            switch block
                % 1st baseline
                case 1
                    stimChans = [5 6];
                    
                    
                    % 1st baseline
                case 2
                    stimChans = [4 3];
                    
                    % 1st post stim
                case 3
                    stimChans = [5 6];
                    
                    % 1st post stim
                case 4
                    stimChans = [4 3];
                    
                    % 2nd baseline
                case 5
                    stimChans = [5 6];
                    
                    % 2nd baseline
                case 6
                    stimChans = [4 3];
                    
                    % 2nd post stim
                case 7
                    stimChans = [5 6];
                    
                    % 2nd post stim
                case 8
                    % they were stimulating the entire time
                    stimChans = [4 3];
                    
                    % 2nd post stim
                case 9
                    % this was a second "post measurement"
                    stimChans = [4 3];
                    
                    % third baseline
                case 10
                    % left side DBS connected
                    stimChans = [5 6];
                    
                    % third baseline
                case 11
                    stimChans = [ 4 3];
                    
                    % 3rd post stim
                case 12
                    stimChans = [5 6];
                    
                    % 3rd post stim
                case 13
                    stimChans = [4 3];
                    
            end
            
        case '3d413'
            tBegin = 2.1; % ms was 2.5 for 42
            tEnd = 35; % ms % was 35 for 426
            
            switch block
                % 1st baseline - patient was under/waking up
                case 1
                    stimChans = [7 8];
                    
                    % 1st baseline - patient was under/waking up
                case 2
                    stimChans = [6 5];
                    
                    % 2nd baseline, patient awake
                case 3
                    stimChans = [7 8];
                    
                    % 2nd baseline, patient awake
                    % pain meds added
                case 4
                    stimChans = [6 5];
                    
                    % 3nd baseline, patient awake, during CT scan
                case 5
                    stimChans = [7 8];
                    
                    % 3nd baseline, patient awake
                case 6
                    stimChans = [6 5];
                    
                    %4th baseline, going under
                case 7
                    stimChans = [7 8];
                    
                    %4th baseline, going under
                case 8
                    stimChans = [6 5];
                    
                    %5th baseline, under
                case 9
                    stimChans = [7 8];
                    
                    %5th baseline, under
                case 10
                    stimChans = [6 5];
                    
            end
            
        case 'fe7df'
            tBegin = 2.1;
            tEnd = 15; %
            
            switch block
                % first baseline - pre stim
                case 1
                    stimChans = [7 8];
                    
                    % second baseline - pre stim
                case 2
                    stimChans = [7 8];
                    
                    % third baseline - pre stim
                case 3
                    stimChans = [7 8];
                    
                    % post 25 ms conditioning - fluid warmer weird noise maybe
                case 4
                    stimChans = [7 8];
                    
                    % post 25 ms conditioning , two minutes after above
                case 5
                    stimChans = [7 8];
                    
                    % post 200 ms conditioning
                case 6
                    stimChans = [7 8];
                    
                    % post 200 ms conditioning
                case 7
                    stimChans = [7 8];
                    
                    % post 200 ms conditioning
                case 8
                    stimChans = [7 8];
                    
                    % post 200 ms conditioning
                case 9
                    stimChans = [7 8];
                    
            end
            
        case 'e6f3c'
            tBegin = 2.1;
            tEnd = 50; %
            
            switch block
                % first baseline - hitting rails
                case 1
                    stimChans = [8 7];
                    
                    % first baseline - hitting rails
                case 2
                    stimChans = [8 7];
                    
                    % first baseline - pre stim - 1.5 2 2.5 0.3 mA whoops
                case 3
                    stimChans = [6 5];
                    
                    % first baseline - pre stim - 1.5 2 2.5 0.3 mA whoops
                case 4
                    stimChans = [6 5];
                    
                    % first baseline - look at this for 6/5 - first one
                    % 1.5 2 2.5 3 mA
                case 5
                    stimChans = [6 5];
                    
                    % first baseline - look at this for 8/7 - first one
                case 6
                    stimChans = [8 7];
                    
                    % first baseline - second one
                case 7
                    stimChans = [8 7];
                    
                    % first baseline - second one
                case 8
                    stimChans = [6 5];
                    
                    % post 200 ms conditioning
                case 9
                    stimChans = [8 7];
                    
                    % post 200 ms conditioning
                case 10
                    stimChans = [6 5];
                    
                    
                    % second baseline
                case 11
                    stimChans = [8 7];
                    
                    % second baseline
                case 12
                    stimChans = [6 5];
                    
                    % post A/A 25 ms conditioning
                case 13
                    stimChans = [8 7];
                    
                    % post A/A 25 ms conditioning
                case 14
                    stimChans = [6 5];
                    
                    % post A/B 25 ms conditioning
                case 15
                    stimChans = [8 7];
                    
                    % post A/B 25 ms conditioning -noisy
                case 16
                    stimChans = [6 5];
                    
                    % post A/B 25 ms conditioning - less noise we hope
                case 17
                    stimChans = [6 5];
                    
            end
            
        case '9f852'
            tBegin = 3;
            tEnd = 25;
            switch block
                % baseline pre stim 1
                case 1
                    stimChans = [5 6];
                    
                    % baseline pre stim 2
                case 2
                    stimChans = [5 6];
                    
                    % baseline post 25 ms A/B - 1
                case 3
                    stimChans = [5 6];
                    
                    % baseline post 25 ms A/B - 2
                case 4
                    stimChans = [5 6];
                    
                    % baseline post 25 ms A/A - 1
                case 5
                    stimChans = [5 6];
                    
                    % baseline post 25 ms A/A - 2
                case 6
                    stimChans = [5 6];
                    
                    % baseline post 200 ms A/B - 1
                case 7
                    stimChans = [5 6];
                    
                    % baseline post 200 ms A/B - 2
                case 8
                    stimChans = [5 6];
                    
                    % baseline post 200 ms A/B - 3
                case 9
                    stimChans = [5 6];
                    % baseline post 200 ms A/B - 4
                case 10
                    stimChans = [5 6];
                    
                    % baseline post 25 ms A/B - 1 - second time
                case 11
                    stimChans = [5 6];
                    
                    % baseline post 25 ms A/B - 2 - second time
                case 12
                    stimChans = [5 6];
                    
            end
            
            
            
        otherwise
            error('unknown sid')
    end
    load(fullfile(SUB_DIR,sid,matlab_dir,experiment,['EP_Measure-' num2str(block) '.mat']));
    
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
    postSamps = round(1000*ECoGfs/1e3);
    tEpoch = [-preSamps:postSamps-1]/ECoGfs*1e3;
    
    goodVec = logical(ones(size(ECoG,2),1));
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
    
    
    % get peak to peak values
    rerefMode = 'none';
    smooth = 1;
    
    [signalPP,pkLocs,trLocs] =  extract_PP_peak_to_peak(stimLevelUniq,epochsEP,tEpoch,stimChans,tBegin,tEnd,rerefMode,[],smooth);
    
    signalPPblock{blockCount} = signalPP;
    pkLocsBlock{blockCount} = pkLocs;
    trLocsBlock{blockCount} = trLocs;
    
    [signalPP,pkLocs,trLocs] =  extract_PP_peak_to_peak_single_trial(stimLevelUniq,epochsEP,tEpoch,stimChans,tBegin,tEnd,rerefMode,[],smooth);
    
    signalPPblockST{blockCount} = signalPP;
    pkLocsBlockST{blockCount} = pkLocs;
    trLocsBlockST{blockCount} = trLocs;
    
    blockCount = blockCount + 1;
    
    clearvars badTrials badTrialLocations
    
    fprintf(['finished block '  num2str(block) '\n'])
    
end


