%% iterate through blocks

blockCount = 1; % do not change, counter variable
blockLabel = {};

if plotCondAvg
    
    figStimBox = figure;
    figCurrent = figure;
    
end
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
            tBegin = 2;
            tEnd = 30;
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
            
            tBegin = 1; % ms
            tEnd = 25; % ms
            
            switch block
                % pre first time
                case 1
                    stimChans = [6 5];
                case 2
                    stimChans = [7 8];
                    badTrials = 1;
                    badTrialLocations = [1:16, 39, 40, 83];
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
            
        case '9eec7'
            
            error('no measurements for this subject')
            
            
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
            
            tBegin = 7;
            tEnd = 70;
            switch block
                % baseline pre stim 1
                case 1
                    stimChans = [5 6];
                    
                    % baseline pre stim 2
                case 2
                    stimChans = [5 6];
                    
                    %  post 25 ms A/B - 1
                case 3
                    stimChans = [5 6];
                    
                    % baseline post 25 ms A/B - 2
                case 4
                    stimChans = [5 6];
                    
                    %  post 25 ms A/A - 1
                case 5
                    stimChans = [5 6];
                    
                    % baseline post 25 ms A/A - 2
                case 6
                    stimChans = [5 6];
                    
                    %post 200 ms A/B - 1
                case 7
                    stimChans = [5 6];
                    
                    % baseline post 200 ms A/B - 2
                case 8
                    stimChans = [5 6];
                    
                    % post 200 ms A/B - 3
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
            
        case '71c6c'
            error('no measurements for this subject')
            
        case '8e907'
            tBegin = 1.8; % ms was 2.5 for 42
            tEnd = 17; % ms % was 35 for 426
            
            switch block
                % baseline pre stim 1
                
                case 1
                    stimChans = [6 7];
                    
                    % baseline pre stim 2, much more drilling
                case 2
                    stimChans = [6 7];
                    
                    % baseline post 200 ms A/B, 500 us pulse width ignore
                case 3
                    stimChans = [6 7];
                    
                    % baseline post 200 ms A/B, 500 us pulse width ignore
                case 4
                    stimChans = [6 7];
                    
                    % baseline post 200 ms A/B, 200 us pulse width
                case 5
                    stimChans = [6 7];
                    
                    % baseline post 200 ms A/B, 200 us pulse width , after
                    % waiting two minutes
                case 6
                    stimChans = [6 7];
                    
                    % baseline post 200 ms A/B, 200 us pulse width, final
                    % one
                case 7
                    stimChans = [6 7];
                    
            end
            
        case '08b13'
            
            tBegin = 3.5; % ms
            tEnd =28; % ms
            switch block
                % baseline pre stim 1
                case 1
                    stimChans = [7 8];
                    
                    % baseline pre stim 1
                case 2
                    stimChans = [5 6];
                    
                    % baseline pre stim 2
                case 3
                    stimChans = [7 8];
                    
                    % baseline pre stim 2 - patient starting to go to sleep
                case 4
                    stimChans = [5 6];
                    
                    % baseline post 200 ms A/B 200 ms delay
                case 5
                    stimChans = [7 8];
                    
                    % second baseline post 200 ms A/B , second baseline to compare
                    % against for the following A/A tests
                    
                case 6
                    stimChans = [7 8];
                    
                    % post 200 ms A/A
                case 7
                    stimChans = [7 8];
                    
                    % post 200 ms A/A - second one
                case 8
                    stimChans = [7 8];
            end
            
        case 'e9c9b'
            
            tBegin = 3; % ms
            tEnd = 55; % ms
            switch block
                % baseline pre stim 1
                case 1
                    stimChans = [7 8];
                    
                    % baseline pre stim 2
                case 2
                    stimChans = [7 8];
                    
                    % post A/A 200
                case 3
                    stimChans = [7 8];
                    
                    % baseline 3 post A/A 200
                case 4
                    stimChans = [7 8];
                    
                    % post A/B 200 ms
                case 5
                    stimChans = [7 8];
                    
                    % baseline post A/B 200
                case 6
                    stimChans = [7 8];
                    
                    % post A/B 25 ms
                case 7
                    stimChans = [7 8];
                    
                    % baseline 5 post A/B
                case 8
                    stimChans = [7 8];
                    
                    %baseline again - 6
                case 9
                    stimChans = [7 8];
                    % post A/A 25 ms
                case 10
                    stimChans = [7 8];
            end
            
        case '41a73'
            
            tBegin = 2.5; % ms
            tEnd = 35; % ms
            switch block
                % baseline pre stim 1
                case 1
                    stimChans = [6 7];
                    
                    % baseline pre stim 2
                case 2
                    stimChans = [6 7];
                    
                    % post A/B 200 ms
                case 3
                    stimChans = [6 7];
                    badTrials = 1;
                    badTrialLocations = [103 196 201:222];
                    
                    % baseline post A/B 200 - baseline 3
                case 4
                    stimChans = [6 7];
                    badTrials = 1;
                    badTrialLocations = [80 132 145:182];
                    %
                    % baseline 4
                case 5
                    stimChans = [6 7];
                    
                    % post A/A 200 ms
                case 6
                    stimChans = [6 7];
                    badTrials = 1;
                    badTrialLocations = [193:240];
                    
                    % baseline post A/A 200 ms - baseline 5
                case 7
                    stimChans = [6 7];
                    
                    % baseline 6
                case 8
                    stimChans = [6 7];
                    
                    %baseline 7
                case 9
                    stimChans = [6 7];
                    
                    % baseline 8
                case 10
                    stimChans = [6 7];
                    
                    % baseline 9
                case 11
                    stimChans = [6 7];
                    
                    % baseline 10
                case 12
                    stimChans = [6 7];
            end
            
        case '68574'
            
            tBegin = 2.5; % ms
            tEnd = 60; % ms
            switch block
                % baseline pre stim 1
                case 1
                    stimChans = [5 6];
                    
                    % baseline pre stim 2
                case 2
                    stimChans = [5 6];
                    
                    % post A/A 100 ms
                case 3
                    stimChans = [5 6];
                    badTrials = 1;
                    badTrialLocations = [74:167];
                    
                    % baseline post A/A 100 ms
                case 4
                    stimChans = [5 6];
                    
                    % post A/B 100
                case 5
                    stimChans = [5 6];
                    
                    % baseline4-  post A/B
                case 6
                    stimChans = [5 6];
                    badTrials = 1;
                    badTrialLocations = [108:156];
                    
                    
                    % post A/A 200 ms
                case 7
                    stimChans = [5 6];
                    badTrials = 1;
                    badTrialLocations = [1:9 43:56];
                    
                    % baseline 5 (post A/A )
                case 8
                    stimChans = [5 6];
                    
                    % post A/B 200 ms
                case 9
                    stimChans = [5 6];
                    
                    % baseline 6 - post A/B 200 ms
                case 10
                    stimChans = [5 6];
                    
                    % baseline 7 - pre DBS
                    
                case 11
                    stimChans = [5 6];
                    badTrials = 1;
                    badTrialLocations = [62 97];
                    
                    % during DBS
                case 13
                    stimChans = [5 6];
                    
                    % post DBS
                case 14
                    stimChans = [5 6];
                    
                    % post DBS 2
                case 15
                    stimChans = [5 6];
            end
            
        case '01fee'
            
            tBegin = 2; % ms
            tEnd = 65; % ms
            switch block
                % baseline pre stim 1
                case 1
                    stimChans = [6 7];
                    
                    % baseline pre stim 2
                case 2
                    stimChans = [6 7];
                    
                    % post A/B 100 ms
                case 3
                    stimChans = [6 7];
                    badTrials = 1;
                    badTrialLocations = [107:129 237:240];
                    % baseline post A/B 100 ms
                case 4
                    stimChans = [6 7];
                    badTrials = 1;
                    badTrialLocations = [107:129 168:240];
                    % post A/A 100
                case 5
                    stimChans = [6 7];
                    
                    % baseline 4-  post A/A
                case 6
                    stimChans = [6 7];
                    
                    % post A/B 200 ms
                case 7
                    stimChans = [6 7];
                    
                    % baseline 5 (post A/B )
                case 8
                    stimChans = [6 7];
                    badTrials = 1;
                    badTrialLocations = [64 88 167];
                    
                    %  during DBS
                case 9
                    stimChans = [6 7];
                    
                    % after DBS 1
                case 10
                    stimChans = [6 7];
                    badTrials = 1;
                    badTrialLocations = [169 204 205];
                    % after DBS 2
                case 11
                    stimChans = [6 7];
                    badTrials = 1;
                    badTrialLocations = [132];
            end
            
        case 'a23ed'
            
            tBegin = 3; % ms
            tEnd = 30; % ms
            switch block
                % baseline pre stim 1
                case 1
                    stimChans = [6 7];
                    badTrials = 1;
                    badTrialLocations = [195:219 239:240];
                    
                    % baseline pre stim 2
                case 2
                    stimChans = [6 7];
                    badTrials = 1;
                    badTrialLocations = [128:161 216:240];
                    
                    % post A/B 200 ms - only 5 minutes!
                case 3
                    stimChans = [6 7];
                    badTrials = 1;
                    badTrialLocations = [62 180 216:240];
                    
                    % baseline 3 post A/B 200 ms  - 15 mins
                case 4
                    stimChans = [6 7];
                    badTrials = 1;
                    badTrialLocations = [78 190];
                    
                    % baseline 4
                case 5
                    stimChans = [6 7];
                    badTrials = 1;
                    badTrialLocations = [221];
                    
                    % post A/B 200 ms - 15 minutes
                case 6
                    stimChans = [6 7];
                    
                    % baseline 5
                case 7
                    stimChans = [6 7];
                    
                    % post A/A 200 ms - 15 minutes - DSB was on during
                    % conditioning
                case 8
                    stimChans = [6 7];
                    
                    % baseline 6 - pre DBS
                case 9
                    stimChans = [6 7];
                    badTrials = 1;
                    badTrialLocations = [78];
            end
            
        otherwise
            error('unknown sid')
    end
    load(fullfile(SUB_DIR,sid,matlab_dir,experiment,['EP_Measure-' num2str(block) '.mat']));
    
    % default values for t begin
    if ~exist('tBegin','var')
        tBegin = 1.8; % ms was 2.5 for 42
        tEnd = 35; % ms % was 35 for 426
    end
    
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
    
    if screenBadChans
        figure
        index = 1;
        stimCommandTimesNew = round((29+stimCommandTimes)*fac);
        epochTemp = getEpochSignal(ECoG,stimCommandTimesNew-preSamps,stimCommandTimesNew+postSamps);
        for jj=1:size(epochTemp,3)
            subplot(15,16,index)
            plot(tEpoch,epochTemp(:,5,jj))
            ylim([-1e-3 1e-3])
            title(num2str(index))
            index = index + 1;
        end
        suptitle(['Block ' num2str(block)])
    end
    
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
    for ii = stimLevelUniqAll
        if sum(stimCommandTimes(stimLevelCommandTimes == ii)) == 0
            stimLevelUniq(count) = nan;
        end
        count = count + 1;
    end
    
    stimLevelUniq = stimLevelUniq(~isnan(stimLevelUniq));
    
    preSamps = round(100*ECoGfs/1e3);
    postSamps = round(490*ECoGfs/1e3);
    tEpoch = [-preSamps:postSamps-1]/ECoGfs*1e3;
    
    % seems to be 29 sample delay between stim command and when the ECoG
    % data starts to move
    % and measured peak
    count = 1;
    
    for ii = stimLevelUniq
        tempStimTimes = round((29+stimCommandTimes(stimLevelCommandTimes == ii))*fac);
        
        % make sure don't exceed edges of trial
        tempEnd = tempStimTimes + postSamps;
        badInds = tempEnd>size(ECoG,1);
        tempStimTimes(badInds) = [];
        stimLevelCell{count} = tempStimTimes;
        count = count + 1;
        
    end
    
    %%
    epochsEP = {};
    
    
    goodVec = logical(ones(size(ECoG,2),1));
    goodVec(stimChans) = 0;
    chansList = [1:8];
    chans = chansList(goodVec);
    ECoG(:,stimChans,:) = 0;
    
    count = 1;
    
    chanInt = chanIntList(1);
    
    for ii = stimLevelUniq
        
        epochUnNormed = getEpochSignal(ECoG,stimLevelCell{count}-preSamps,stimLevelCell{count}+postSamps);
        epochNormed =  epochUnNormed-repmat(mean(epochUnNormed(tEpoch<0,:,:),1), [size(epochUnNormed, 1), 1]);
        epochsEP{count} = epochNormed;
        if plotCondAvg
            smallMultiples_DBS(mean(epochsEP{count},3),tEpoch/1e3,'type2',stimChans);
            title(['stim pair = '  num2str(stimChans) ' Stim Level ' num2str(ii)])
            
            if savePlot
                SaveFig(OUTPUT_DIR, [sid '_EP_block_' num2str(block) '_stimChans_' num2str(stimChans(1)) '_' num2str(stimChans(2)) '_stimLev_' num2str(ii)], 'png', '-r600');
            end
        end
        
        %%
        [processedSig,templateDictCell,templateTrial,startInds,endInds] = analyFunc.template_subtract(epochUnNormed,'type',type,...
            'fs',ECoGfs,'plotIt',plotIt,'pre',pre,'post',post,'stimChans',stimChans,...
            'useFixedEnd',useFixedEnd,'fixedDistance',fixedDistance,...
            'distanceMetricDbscan',distanceMetricDbscan,'distanceMetricSigMatch',distanceMetricSigMatch,...
            'recoverExp',recoverExp,'normalize',normalize,'amntPreAverage',amntPreAverage,...
            'minDuration',minDuration,'bracketRange',bracketRange,'threshVoltageCut',threshVoltageCut,...
            'threshDiffCut',threshDiffCut,'onsetThreshold',onsetThreshold,'chanInt',chanInt,...
            'minPts',minPts,'minClustSize',minClustSize,'outlierThresh',outlierThresh,'useProcrustes',useProcrustes);
        %%
        % visualization
        % of note - more visualizations are created here, including what the
        % templates look like on each channel, and what the discovered templates are
        xlims = [-50 400];
        %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        average = 1;
        modePlot = 'avg';
        ylims = [-0.6 0.6];
        vizFunc.small_multiples_time_series(processedSig,tEpoch,'type1',stimChans,'type2',0,'xlims',xlims,'ylims',ylims,'modePlot',modePlot,'highlightRange',trainDuration)
        
        %%%%%% wavelet
        fprintf(['-------Beginning wavelet analysis-------- \n'])
        
        timeRes = 0.01; % 10 ms bins
        
        [powerout,fMorlet,tMorlet,~] = analyFunc.waveletWrapper(processedSig,ECoGfs,timeRes,stimChans);
        %
        fprintf(['-------Ending wavelet analysis-------- \n'])
        
        % additional parameters
        postStim = max(tEpoch);
        sampsPostStim = round(postStim/1e3*ECoGfs);
        
        preStim = min(tEpoch);
        sampsPreStim = round(preStim/1e3*ECoGfs);
        
        tMorlet = linspace(preStim,postStim,length(tMorlet))/1e3;
        % normalize data
        dataRef = powerout(:,tMorlet<-0.005 & tMorlet>-0.08,:,:);
        dataRef = powerout(:,:,:,:);

        %
        [normalizedData] = analyFunc.normalize_spectrogram_wavelet(dataRef,powerout);
        individual = 1;
        average = 1;
        
        for chanInt = chanIntList
            vizFunc.visualize_wavelet_channel_no_raw(normalizedData,tMorlet,fMorlet,processedSig,...
                tEpoch/1e3,chanInt,individual,average,xlimsWave)
        end
        
        for chanInt = chanIntList
            vizFunc.visualize_wavelet_channel(normalizedData,tMorlet,fMorlet,processedSig,...
                tEpoch/1e3,epochUnNormed,chanInt,individual,average,xlimsWave)
        end
        
        for chanInt = chanIntList
            vizFunc.visualize_wavelet_channel_no_raw_not_normalized(powerout,tMorlet,fMorlet,processedSig,...
                tEpoch/1e3,chanInt,individual,average,xlimsWave)
        end
        %         %
        %         for chanInt = chanIntList
        %             vizFunc.visualize_wavelet_channel_small(normalizedData,tMorlet,fMorlet,processedSig,...
        %                 tEpoch,dataInt,chanInt,individual,average)
        %       end
        %
        ylimsSpect = [5 300];
        vizFunc.small_multiples_spectrogram(normalizedData,tMorlet,fMorlet,'type1',stimChans,'type2',0,'xlims',xlims,'ylims',ylimsSpect);
        
        
        count = count + 1;
        
    end
    
    %%
    epochsEPblock{blockCount} = epochsEP;
    
    % get peak to peak values
    smooth = 1;
    
    % example channel for peakt o peak
    chanIntTemp = chanIntList(1);
    
    
    [signalPPAvg,pkLocsAvg,trLocsAvg] =  extract_PP_peak_to_peak(stimLevelUniq,epochsEP,tEpoch,stimChans,...
        tBegin,tEnd,rerefMode,chanReref,smooth,chanIntTemp);
    
    signalPPblock{blockCount} = signalPPAvg;
    pkLocsBlock{blockCount} = pkLocsAvg;
    trLocsBlock{blockCount} = trLocsAvg;
    
    
    [signalPP,pkLocs,trLocs] =  extract_PP_peak_to_peak_single_trial(stimLevelUniq,epochsEP,tEpoch,...
        stimChans,tBegin,tEnd,rerefMode,chanReref,smooth,avgTrials,numAvg,chanIntTemp);
    
    signalPPblockST{blockCount} = signalPP;
    pkLocsBlockST{blockCount} = pkLocs;
    trLocsBlockST{blockCount} = trLocs;
    
    [signalPPfromAvg,pkLocsFromAvg,trLocsFromAvg] =  extract_PP_ind_from_avg_max_min(stimLevelUniq,epochsEP,...
        tEpoch,tBegin,tEnd,pkLocsAvg,trLocsAvg,stimChans,smooth,rerefMode,chanReref,avgTrials,numAvg);
    
    signalPPblockSTfromAvg{blockCount} = signalPPfromAvg;
    pkLocsBlockSTfromAvg{blockCount} = pkLocsFromAvg;
    trLocsBlockSTfromAvg{blockCount} = trLocsFromAvg;
    
    % make labels to keep track of each trial
    for counter = 1:size(signalPP,2)
        blockLabel{blockCount}{counter} = repmat(stimLevelUniq(counter),size(signalPP{counter},2),1);
    end
    
    
    blockCount = blockCount + 1;
    
    clearvars badTrials badTrialLocations
    
    %     figure(figStimBox)
    %     hold on
    %     plot(stimBox(:,1))
    %
    %     figure(figCurrent)
    %     hold on
    %     plot(stimProgrammed(:,1))
    %
    %     fprintf(['finished block '  num2str(block) '\n'])
    %
end
% figure(figStimBox)
%
% legend
%
% figure(figCurrent)
% legend


