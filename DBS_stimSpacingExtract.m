%% DJC - 8-29-2016 - DBS analysis script, updated 5/23/2017
% This is to extract the neural data

%% initialize output and meta dir
% clear workspace - be in the directory with all scripts necessary
close all; clear all; clc
% add path for scripts to work with data tanks
addpath('.\scripts')

% set path, set
Z_ConstantsDBS

% subject directory, change as needed
% for David
SUB_DIR = fullfile(myGetenv('subject_dir'));

%% load in subject

% this is from my z_constants

% for param sweep, look at subjects 1,2,9
%sid = input('what is the sid?\n','s');
%sid = SIDS{2}; % MUST SWITCH THIS, either 1,2,9
sid = '50ad9'
% load in tank
switch sid
    case 'bb908'
        
        [structureData,filepath] = promptForTDTrecording;
        split_path = split(filepath,"\");
        fileName = split_path{end};
        
        Sing = structureData.Sing;
        Stim = structureData.Stim;
        Valu = structureData.Valu;
        Cond = structureData.Cond;
        DBSs = structureData.DBSs;
        ECOG = structureData.ECOG;
        
        dbsElectrodes = DBSs.data;
        dbs_fs = DBSs.info.SamplingRateHz;
        
        ECOGelectrodes = ECOG.data;
        ECOG_fs = ECOG.info.SamplingRateHz;
        
        stimBox = Stim.data;
        stim_fs = Stim.info.SamplingRateHz;
        
        stimProgrammed = Sing.data;
        
        stimSampDeliver = Cond.data(:,1);
        condition = Cond.data(:,2);
        ttlPulse = Cond.data(:,3);
        cond_fs = Cond.info.SamplingRateHz;
        
        
        
        
    case '80301'
        [structureData,filepath] = promptForTDTrecording;
        split_path = split(filepath,"\");
        fileName = split_path{end};
        Sing = structureData.Sing;
        Stim = structureData.Stim;
        Valu = structureData.Valu;
        Cond = structureData.Cond;
        DBSs = structureData.DBSs;
        ECOG = structureData.ECOG;
        
        dbsElectrodes = DBSs.data;
        dbs_fs = DBSs.info.SamplingRateHz;
        
        ECOGelectrodes = ECOG.data;
        ECOG_fs = ECOG.info.SamplingRateHz;
        
        stimBox = Stim.data;
        stim_fs = Stim.info.SamplingRateHz;
        
        stimProgrammed = Sing.data;
        
        stimSampDeliver = Cond.data(:,1);
        condition = Cond.data(:,2);
        ttlPulse = Cond.data(:,3);
        cond_fs = Cond.info.SamplingRateHz;
        
        % deal with that chunk in stimParam12 that's bad
        if strcmp(fileName,'paramsweep-12.mat')
            BadData = ones(size(dbsElectrodes,1),1);
            BadData(5.983e6:6.84e6) = 0;
            BadData = logical(BadData);
            dbsElectrodes = dbsElectrodes(BadData,:);
            ECOGelectrodes = ECOGelectrodes(BadData,:);
            stimBox = stimBox(BadData,:);
            stimProgrammed = stimProgrammed(BadData,:);
            stimSampDeliver = stimSampDeliver(BadData,:);
            condition = condition(BadData,:);
            ttlPulse = ttlPulse(BadData,:);
        end
        
    case '1dd75'
        [structureData,filepath] = promptForTDTrecording;
        split_path = split(filepath,"\");
        fileName = split_path{end};
        Sing = structureData.Sing;
        Stim = structureData.Stim;
        Valu = structureData.Valu;
        Cond = structureData.Cond;
        DBSs = structureData.DBSs;
        ECOG = structureData.ECOG;
        
        dbsElectrodes = DBSs.data;
        dbsElectrodes = dbsElectrodes(:,1:4); % only recorded left side DBS
        dbs_fs = DBSs.info.SamplingRateHz;
        
        % 1->8 are ECoG, 9->12 are cerebellar
        ECOGelectrodes = ECOG.data;
        ECOG_fs = ECOG.info.SamplingRateHz;
        
        stimBox = Stim.data;
        stim_fs = Stim.info.SamplingRateHz;
        
        stimProgrammed = Sing.data;
        
        stimSampDeliver = Cond.data(:,1);
        %condition = Cond.data(:,2); % this condition did not save right
        %for this subject! Will work from now on
        
        % from DBS_10_20_2016_condition_1.txt, have to sync up with stim
        % samp deliver
        condition_file = [2,3,4,1,4,3,1,2,1,4,2,3,2,4,1,3,1,2,4,3,1,4,3,2,2,4,3,1,4,3,1,2,4,1,2,3,3,4,2,1,1,3,4,2,1,3,2,4,4,3,1,2,3,1,4,2,3,4,1,2];
        
        ttlPulse = Cond.data(:,3);
        cond_fs = Cond.info.SamplingRateHz;
        
    case '50ad9'
        
        [structureData,filepath] = promptForTDTrecording;
        split_path = split(filepath,"\");
        fileName = split_path{end};
        
        Sing = structureData.Sing;
        Stim = structureData.Stim;
        Valu = structureData.Valu;
        Cond = structureData.Cond;
        DBSs = structureData.DBSs;
        ECOG = structureData.ECOG;
        
        dbsElectrodes = DBSs.data;
        dbs_fs = DBSs.info.SamplingRateHz;
        
        ECOGelectrodes = ECOG.data;
        ECOG_fs = ECOG.info.SamplingRateHz;
        
        stimBox = Stim.data;
        stim_fs = Stim.info.SamplingRateHz;
        
        stimProgrammed = Sing.data;
        
        stimSampDeliver = Cond.data(:,1);
        condition = Stim.data(:,2); % note - this is changed from the other previous subjects
        ttlPulse = Cond.data(:,3);
        cond_fs = Cond.info.SamplingRateHz;
    case '695e1'
        
        [structureData,filepath] = promptForTDTrecording;
        split_path = split(filepath,"\");
        fileName = split_path{end};
        
        Sing = structureData.Sing;
        Stim = structureData.Stim;
        Valu = structureData.Valu;
        Cond = structureData.Cond;
        DBSs = structureData.DBSs;
        ECOG = structureData.ECOG;
        
        dbsElectrodes = DBSs.data;
        dbs_fs = DBSs.info.SamplingRateHz;
        
        ECOGelectrodes = ECOG.data;
        ECOG_fs = ECOG.info.SamplingRateHz;
        
        stimBox = Stim.data;
        stim_fs = Stim.info.SamplingRateHz;
        
        stimProgrammed = Sing.data;
        
        stimSampDeliver = Cond.data(:,1);
        condition = Stim.data(:,2); % note - this is changed from the other previous subjects
        ttlPulse = Cond.data(:,3);
        cond_fs = Cond.info.SamplingRateHz;
        
end

%% decide what to plot

% ui box for input
prompt = {'Plot stimulation monitor and current to be delivered (time series?) y or n ',...
    'Plot time series of DBS and ECoG electrodes? y or n','Plot Specific channels or conditions of interest? y or n',...
    'Do voltage analysis? y or n',...
    'Find stim delivery peaks & Plot histogram of DBS and ECoG electrodes? y or n',...
    'Find cceps? y or n',...
    'Plot CCEPs','Save output file','Save internal CCEPs','peak verification','Subtract pre period'};
dlg_title = 'Channel of Interest';
num_lines = 1;
defaultans = {'n','y','y','n','n','y','y','n','y','y','n'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

plotStim = answer{1};
plotTime = answer{2};
plotCond = answer{3};
voltageAnalysis = answer{4};
plotHist = answer{5};
ccepAnalysis = answer{6};
plotCCEP = answer{7};
saveOutput = answer{8};
saveOutputInternal = answer{9};
assurePeaks = answer{10};
subtractPre = answer{11};
%% plot stim
%
if strcmp(plotStim,'y')
    figure
    hold on
    for i = 1:size(stimBox,2)
        
        t = (0:length(stimBox)-1)/stim_fs;
        subplot(2,2,i)
        plot(t*1e3,stimBox(:,i))
        title(sprintf('Channel %d',i))
        
        
    end
    
    
    xlabel('Time (ms)')
    ylabel('Amplitude (V)')
    
    %subtitle('Stimulation Channels')
end

%% Sing is wave to be delivered

% build a burst table with the timing of stimuli from the stim file
bursts = [];

Sing1 = stimProgrammed(:,1);


% trying something like A_BuildStimTables from BetaStim


stimMask = stimSampDeliver~=0;


% sample length of train - 500 ms
sampsEnd = round(0.5*stim_fs); % CHANGED THIS TO ROUND RATHER THAN FLOOR? DOESNT MATTER HERE PROBABLY 

bursts(2,:) = find(stimMask==1);
bursts(3,:) = bursts(2,:) + repmat(sampsEnd,size(bursts(2,:)));

stims = squeeze(getEpochSignal(Sing1,(bursts(2,:)-1),(bursts(3,:))+1));
t = (0:size(stims,1)-1)/stim_fs;
t = t*1e3;
if strcmp(plotStim,'y')
    
    figure
    % for first subject, stim_epcohed 1 seems to be off by a sample
    plot(t,stims(:,2:end))
    xlabel('Time (ms)');
    ylabel('Voltage to be delivered')
    title('Voltage to be delivered')
    
    
end
%% Plot stims with info from above

stim1 = stimBox(:,1);
stim1Epoched = squeeze(getEpochSignal(stim1,(bursts(2,:)-1),(bursts(3,:))+1));
t = (0:size(stim1Epoched,1)-1)/stim_fs;
t = t*1e3;

if strcmp(plotStim,'y')
    
    figure
    plot(t,stim1Epoched)
    xlabel('Time (ms');
    ylabel('Voltage (V)');
    title('Finding the delay between current output and stim delivery')
    
    hold on
    
    plot(t,stims)
    
    % get the delay in stim times
end

delay = floor(0.1434*stim_fs/1e3);

if strcmp(plotStim,'y')
    
    % plot the appropriately delayed signal
    figure
    stimTimesBegin = bursts(2,:)-1+delay;
    stimTimesEnd = bursts(3,:)+1+delay;
    stim1Epoched = squeeze(getEpochSignal(stim1,stimTimesBegin,stimTimesEnd));
    t = (0:size(stim1Epoched,1)-1)/stim_fs;
    t = t*1e3;
    figure
    plot(t,stim1Epoched(:,2:end))
    hold on
    plot(t,stims(:,2:end))
    xlabel('Time (ms');
    ylabel('Voltage (V)');
    title('Stim voltage monitoring with delay added in')
    
end

%% extract data - NO NEED FOR THIS SINCE DATA IS SAMPLED HIGH HERE
%
% try and account for delay for the stim times
stimTimes = bursts(2,:)-1+delay;


prompt = {'time to look before stimulation (seconds) (If wanting to do stimulation pulse analysis, set to 0) (If wanting to look at internal CCEPs, set to 0.1) (If wanting to do external to train CCEP, -0.495)'...
    ,'Time to look after stimulation signal (seconds) (If wanting to do stimulation pulse analysis, set to 0.495) (If wanting to look at internal CCEPs, set to 0.5) (If wanting to do external to train CCEP, 0.795) '...
    'What was the stimulation frequency?'};
dlg_title = 'How much to analyze';
num_lines = 1;
defaultans = {'0.1','0.5','185'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

pre = str2num(answer{1});
post  = str2num(answer{2});
stimFreq = str2num(answer{3});


% DJC 7-7-2016, changed presamps and post samps to 1 second
presamps = floor(pre * stim_fs); % pre time in sec
postsamps = floor(post * stim_fs); % post time in sec, % modified DJC to look at up to 300 ms after

%% get the data epochs


dataEpochedECOG = squeeze(getEpochSignal(ECOGelectrodes,stimTimes-presamps,stimTimes+postsamps));
dataEpochedDBS = squeeze(getEpochSignal(dbsElectrodes,stimTimes-presamps,stimTimes+postsamps));

numEco = size(dataEpochedECOG,2);
numDBS = size(dataEpochedDBS,2);

% mean subtract

if strcmp(sid,'bb908')
    
    ECoG_ave = mean(dataEpochedECOG,1);
    DBS_ave = mean(dataEpochedDBS,1);
    
    dataEpochedECOG = dataEpochedECOG - repmat(ECoG_ave,size(dataEpochedECOG,1),1,1);
    dataEpochedDBS = dataEpochedDBS - repmat(DBS_ave,size(dataEpochedECOG,1),1,1);
end
%%%%%%%%%%%%%%%%% 80301 wasnt DC coupled
% if strcmp(sid,'80301')
%
%     ECoG_ave = mean(dataEpochedECOG,1);
%     DBS_ave = mean(dataEpochedDBS,1);
%
%     dataEpochedECOG = dataEpochedECOG - repmat(ECoG_ave,size(dataEpochedECOG,1),1,1);
%     dataEpochedDBS = dataEpochedDBS - repmat(DBS_ave,size(dataEpochedECOG,1),1,1);
% end

% set the time vector to be set by the pre and post samps
t = (-presamps:postsamps-1)*1e3/ECOG_fs;

% get conditions
if strcmp(sid,'1dd75')
    ucondition = unique(condition_file);
    %ucondition = ucondition(2:end);
    condition = zeros(length(stimSampDeliver),1);
    condition(stimTimes) = condition_file;
else
    ucondition = unique(condition);
    ucondition = ucondition(2:end);
end

ECoG_sep = {};
DBS_sep = {};

for i = 1:length(ucondition)
    
    ECoG_sep{i} = dataEpochedECOG(:,:,condition(stimTimes)==ucondition(i));
    DBS_sep{i} = dataEpochedDBS(:,:,condition(stimTimes)==ucondition(i));
    
end


%%


prompt = {'Which side was stimulated? L or R ',...
    'Were both DBS leads in? single or both ','1st DBS stim channel (active) '...
    '2nd DBS stim channel (ground) '};
dlg_title = 'Channel of Interest';
num_lines = 1;
defaultans = {'R','both','1','2'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

side = answer{1};
numLeads = answer{2};
stim_chan1 = str2num(answer{3});
stim_chan2 = str2num(answer{4});
stimChans = [stim_chan1 stim_chan2];

if strcmp(saveOutput,'y')
    
    
    save(fullfile(OUTPUT_DIR, ['stim_',side,'_',numLeads,'DBS_' num2str(stim_chan1),'_',num2str(stim_chan2)]),...
        'ucondition','ECOG_fs','dbs_fs','DBS_sep','ECoG_sep','t','presamps','postsamps','pre','post');
    return
    
end

%% plot ECoG Electrodes

if strcmp(plotTime,'y')
    
    for i = 1:length(ucondition)
        
        figure
        ECoG_temp = ECoG_sep{i};
        
        for j = 1:numEco
            
            subplot(4,4,j);
            plot(t,squeeze(ECoG_temp(:,j,:)));
            title(['Channel ' num2str(j)]);
        end
        
        xlabel('time (ms)')
        ylabel('voltage (V)')
        
        subtitle(['ECoG Electrodes, Condition ' num2str(i)]);
        
    end
    
    
    
    %% plot DBS Electrodes
    
    for i = 1:length(ucondition)
        
        figure
        DBS_temp = DBS_sep{i};
        
        for j = 1:numDBS
            
            subplot(4,2,j);
            plot(t,squeeze(DBS_temp(:,j,:)));
            title(['Channel ' num2str(j)]);
        end
        
        xlabel('time (ms)')
        ylabel('voltage (V)')
        
        subtitle(['DBS Electrodes, Condition ' num2str(i)]);
        
    end
    
end
%% plot channel of interest

if strcmp(plotCond,'y')
    % ui box for input
    prompt = {'ECoG Channel of interest?','DBS Channel of interest','Condition of Interest?'};
    dlg_title = 'Channel of Interest';
    num_lines = 1;
    defaultans = {'8','1','4'};
    answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
    
    ecog_chanInt = str2num(answer{1});
    DBS_chanInt = str2num(answer{2});
    cond_int = str2num(answer{3});
    
    ECoG_temp = ECoG_sep{cond_int};
    DBS_temp = DBS_sep{cond_int};
    
    figure
    plot(t,squeeze(ECoG_temp(:,ecog_chanInt,:)))
    
    xlabel('time (ms)')
    ylabel('voltage (V)')
    title(['ECoG Channel ', num2str(ecog_chanInt), ' for Condition ', num2str(cond_int)]);
    
    
    figure
    plot(t,squeeze(DBS_temp(:,DBS_chanInt,:)))
    
    
    xlabel('time (ms)')
    ylabel('voltage (V)')
    title(['DBS Channel ', num2str(DBS_chanInt), ' for Condition ', num2str(cond_int)]);
    
    
    %% look at averages for condition of interest
    DBS_aveCond  = squeeze(mean(DBS_temp,3));
    ECoG_aveCond = squeeze(mean(ECoG_temp,3));
    
    figure
    plot(t,DBS_aveCond);
    xlabel('time (ms)')
    ylabel('voltage (V)')
    title(['Average DBS recording across channels for condition ', num2str(cond_int)])
    
    figure
    plot(t,ECoG_aveCond);
    xlabel('time (ms)')
    ylabel('voltage (V)')
    title(['Average ECoG recording across channels for condition ', num2str(cond_int)])
    
    %% 11-2-2016 - look at subplots of condition of interest - useful for CCEPs after stimulation train ends
    
    
    figure
    for j = 1:numEco
        subplot(4,4,j)
        plot(t,ECoG_aveCond(:,j))
        xlabel('time (ms)')
        ylabel('Voltage (V)')
        title(['Channel ',num2str(j)])
        
    end
    subtitle(['ECoG EP response outside train for condition = ' num2str(cond_int)])
    
    figure
    for j = 1:numDBS
        subplot(2,4,j)
        plot(t,DBS_aveCond(:,j))
        xlabel('time (ms)')
        ylabel('Voltage (V)')
        
        % put a box around the stimulation channels of interest if need be
        if ismember(j,stimChans)
            ax = gca;
            ax.Box = 'on';
            ax.XColor = 'red';
            ax.YColor = 'red';
            ax.LineWidth = 2;
            title(['Channel ',num2str(j)],'color','red');
            
        else
            title(['Channel ',num2str(j)]);
        end
        
    end
    
    subtitle(['DBS EP responses outside train  for condition = ' num2str(cond_int)])
    
    
end

%% get average peaks of waveform
% right now this is particular for looking at the DBS voltage peaks
% delivered

if strcmp(ccepAnalysis,'y') || strcmp(voltageAnalysis,'y')
    
    DBS_peak_pos = {};
    ECoG_peak_pos = {};
    
    DBS_peak_neg = {};
    ECoG_peak_neg = {};
    
    factor = -1;
    ECoG_neg = cellfun(@(x) x*factor,ECoG_sep,'un',0);
    DBS_neg = cellfun(@(x) x*factor,DBS_sep,'un',0);
    
    % ccep
    
    
    if stimFreq == 185
        numPeaks = 92; % empircally found - if 185 Hz
        
    elseif stimFreq == 20
        numPeaks = 9;
    end
    numConds = 15;
    
    numTotal = numPeaks*numConds;
    
    
    preCCEP = floor(3e-3 * stim_fs); % pre time in sec
    
    if stimFreq == 185
        
        postCCEP = floor(6e-3 * stim_fs); % post time in sec
        
    elseif stimFreq == 20
        
        postCCEP = floor(51e-3*stim_fs);
    end
    
    ECoG_sepCCEPinternal = {};
    DBS_sepCCEPinternal = {};
    
    
    % find locations of stimulations using condiiton 4, pick channel of
    % interest
    
    % dbs_7,8,10,11 (:,8,:) is good
    % dbs_9,12 (:,5,:)
    
    prompt = {'Which DBS electrode to use for extract CCEPs?'};
    dlg_title = 'Electrode to extract (for bb908 use 8 for dbs_7,8,10,11 or 5/6 for dbs_9,12 ';
    num_lines = 1;
    defaultans = {'8'};
    answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
    
    chanExtract = str2num(answer{1});
    
    dbs_condP = DBS_sep{cond_int};
    dbs_stackP = squeeze(dbs_condP(:,chanExtract,:));
    dbs_stackP = dbs_stackP(:);
    %%%%%%%%%%%%%%%%%%%% have to tune peak height to extract
    % 11/1/2016 - set min peak height to 2.5e-3 , was 5e-3 before, try
    % 5e-4? -6-6-2017
    % MinPeakH = 2.5e-3;
    MinPeakH = 0.6 * max(dbs_stackP);
    if stimFreq == 185
        [DBS_peakFind_pos,locs] = findpeaks(dbs_stackP,dbs_fs,'MinPeakDistance',0.0048,'NPeaks',numTotal,'MinPeakHeight',MinPeakH);
        findpeaks(dbs_stackP,dbs_fs,'MinPeakDistance',0.0048,'NPeaks',numTotal,'MinPeakHeight',MinPeakH)
    elseif stimFreq == 20
        [DBS_peakFind_pos,locs] = findpeaks(dbs_stackP,dbs_fs,'MinPeakDistance',0.049,'NPeaks',numTotal,'MinPeakHeight',MinPeakH);
        findpeaks(dbs_stackP,dbs_fs,'MinPeakDistance',0.049,'NPeaks',numTotal,'MinPeakHeight',MinPeakH)
    end
    locs = round(locs *stim_fs);
    
    ECoG_sepCCEPinternal = {};
    DBS_sepCCEPinternal = {};
    
    % look at CCEPs inside of stimulation window
    
    presampsCCEP =  floor(0.003*dbs_fs);
    if stimFreq == 185
        
        postsampsCCEP = floor(0.005*dbs_fs);
    elseif stimFreq == 20
        postsampsCCEP = floor(0.05*dbs_fs);
        
    end
    % set the time vector to be set by the pre and post samps
    tCCEP = (-presampsCCEP:postsampsCCEP-1)*1e3/ECOG_fs;
    
end

%%
if strcmp(ccepAnalysis,'y') || strcmp(voltageAnalysis,'y')
    
    for i = 1:length(ucondition)

        % differing numbers of trials
        if strcmp(fileName,'paramsweep-12.mat')
            i =4;
        end
        
        dbs_condP = DBS_sep{i};
        dbs_condN = DBS_neg{i};
        ECoG_condP = ECoG_sep{i};
        ECoG_condN = ECoG_neg{i};
        %%%%%%%%%%%%%%%%%%%% have to tune peak height to extract
        
        % trials to exclude? - trial 14 for 50ad9 condition 3 paramsweep-6
        if i == 3 && strcmp('50ad9',sid)
        trialsExclude = [14];
        mask = ones(size(dbs_condP,3),1);
        mask(trialsExclude) = 0;
        mask = logical(mask);
        else
            mask = ones(size(dbs_condP,3),1);
        end
        
        dbs_stackP = squeeze(dbs_condP(:,chanExtract,mask));
        dbs_stackP = dbs_stackP(:);
        MinPeakH_P = 0.8 * max(dbs_stackP);
        
        dbs_stackN = squeeze(dbs_condN(:,chanExtract,mask));
        dbs_stackN = dbs_stackN(:);
        MinPeakH_N = 0.8 * max(dbs_stackN);
        
        if stimFreq == 185
            [DBS_peakFind_pos,locs_P] = findpeaks(dbs_stackP,dbs_fs,'MinPeakDistance',0.0048,'NPeaks',numTotal,'MinPeakHeight',MinPeakH_P);
            findpeaks(dbs_stackP,dbs_fs,'MinPeakDistance',0.0048,'NPeaks',numTotal,'MinPeakHeight',MinPeakH_P)
        elseif stimFreq == 20
            [DBS_peakFind_pos,locs_P] = findpeaks(dbs_stackP,dbs_fs,'MinPeakDistance',0.049,'NPeaks',numTotal,'MinPeakHeight',MinPeakH_P);
            findpeaks(dbs_stackP,dbs_fs,'MinPeakDistance',0.049,'NPeaks',numTotal,'MinPeakHeight',MinPeakH_P)
        end
        
        if stimFreq == 185
            [DBS_peakFind_neg,locs_N] = findpeaks(dbs_stackN,dbs_fs,'MinPeakDistance',0.0048,'NPeaks',numTotal,'MinPeakHeight',MinPeakH_N);
            findpeaks(dbs_stackN,dbs_fs,'MinPeakDistance',0.0048,'NPeaks',numTotal,'MinPeakHeight',MinPeakH_N)
        elseif stimFreq == 20
            [DBS_peakFind_neg,locs_N] = findpeaks(dbs_stackN,dbs_fs,'MinPeakDistance',0.049,'NPeaks',numTotal,'MinPeakHeight',MinPeakH_N);
            findpeaks(dbs_stackN,dbs_fs,'MinPeakDistance',0.049,'NPeaks',numTotal,'MinPeakHeight',MinPeakH_N)
        end
        
        if assurePeaks
        % help make sure align stims properly
        
        locs_P = round(locs_P *stim_fs);
        locs_N = round(locs_N * stim_fs);
        
        rLocs_P = repmat(locs_P,[1 11]) + [-5:5];
        
        rLocs_N = repmat(locs_N,[1 11]) + [-5:5];
        
        for ind_1 = 1:length(locs_P)
            count = 0;
            
            for  ind_2= 1:length(rLocs_N)
                count_temp = sum(any(locs_P(ind_1)==rLocs_N(ind_2,:),2));
                count = count + count_temp;
            end
            if count == 0
                locs_P(ind_1) = 0;
            end
        end
        
       locs_P(locs_P == 0) =[];
        
        for ind_1 = 1:length(locs_N)
            count = 0;
            
            for  ind_2= 1:length(rLocs_P)
                count_temp = sum(any(locs_N(ind_1)==rLocs_P(ind_2,:),2));
                count = count + count_temp;
            end
            if count == 0
                locs_N(ind_1) = 0;
            end
        end
        
               locs_N(locs_N == 0) =[];

        if locs_P<locs_N
            locs = locs_P;
        else
            locs = locs_N;
        end
        
        % once selected, go with those locs
        else
        locs = floor(locs *stim_fs);
        end
        
        for j = 1:size(dbs_condP,2)
            dbs_stackP = squeeze(dbs_condP(:,j,:));
            dbs_stackP = dbs_stackP(:);
            dbs_stackN = squeeze(dbs_condN(:,j,:));
            dbs_stackN = dbs_stackN(:);
            if strcmp(voltageAnalysis,'y')
                post = 10;
                pre = 5;
                if stimFreq == 185
                    %[DBS_peakFind_pos] = findpeaks(dbs_stackP,dbs_fs,'MinPeakDistance',0.005,'NPeaks',numTotal);
                    dbs_temp = squeeze(getEpochSignal(dbs_stackP,locs-pre,locs+post));
                elseif stimFreq == 20
                    %[DBS_peakFind_pos] = findpeaks(dbs_stackP,dbs_fs,'MinPeakDistance',0.05,'NPeaks',numTotal);
                    dbs_temp = squeeze(getEpochSignal(dbs_stackP,locs-pre,locs+post));
                    
                end
                %DBS_peak_pos{i}{j} = DBS_peakFind_pos;
                DBS_peak_pos{i}{j} = max(dbs_temp,[],1);
                % use locations from this to find CCEP peaks
                
                
                if stimFreq == 185
                    
                    %[DBS_peakFind_neg,~] = findpeaks(dbs_stackN,dbs_fs,'MinPeakDistance',0.005);
                    dbs_temp = squeeze(getEpochSignal(dbs_stackN,locs-pre,locs+post));
                    
                elseif stimFreq == 20
                    %[DBS_peakFind_neg,~] = findpeaks(dbs_stackN,dbs_fs,'MinPeakDistance',0.05);
                    dbs_temp = squeeze(getEpochSignal(dbs_stackN,locs-pre,locs+post));
                    
                    
                end
                %DBS_peak_neg{i}{j} = DBS_peakFind_neg;
                DBS_peak_neg{i}{j} = max(dbs_temp,[],1);
                
            end
            if strcmp(ccepAnalysis,'y')
                dbs_temp = squeeze(getEpochSignal(dbs_stackP,locs-presampsCCEP,locs+postsampsCCEP));
                
                dbs_epoch_ave = mean(dbs_temp((tCCEP<-0.25),:),1);
                                if subtractPre
                dbs_temp = dbs_temp - repmat(dbs_epoch_ave,size(dbs_temp,1),1,1);
                                end
                DBS_sepCCEPinternal{i}(:,:,j) = dbs_temp;
            end
            
        end
        
        for j = 1:size(ECoG_condP,2)
            
            ECoG_stackP = squeeze(ECoG_condP(:,j,:));
            ECoG_stackP = ECoG_stackP(:);
            
            ECoG_stackN = squeeze(ECoG_condN(:,j,:));
            ECoG_stackN = ECoG_stackN(:);
            
            if strcmp(voltageAnalysis,'y')
                post = 10;
                pre = 5;
                
                if stimFreq == 185
                    %[ECoG_peakFind_pos,~] = findpeaks(ECoG_stackP,dbs_fs,'MinPeakDistance',0.005);
                    eco_temp = squeeze(getEpochSignal(ECoG_stackP,locs-pre,locs+post));
                    
                elseif stimFreq == 20
                    %[ECoG_peakFind_pos,~] = findpeaks(ECoG_stackP,dbs_fs,'MinPeakDistance',0.05);
                    eco_temp = squeeze(getEpochSignal(ECoG_stackP,locs-pre,locs+post));
                    
                end
                
                %ECoG_peak_pos{i}{j} = ECoG_peakFind_pos;
                ECoG_peak_pos{i}{j} = max(eco_temp,[],1);
                
                
                if stimFreq == 185
                    %[ECoG_peakFind_neg,~] = findpeaks(ECoG_stackN,dbs_fs,'MinPeakDistance',0.005);
                    eco_temp = squeeze(getEpochSignal(ECoG_stackN,locs-pre,locs+post));
                    
                elseif stimFreq == 20
                    %[ECoG_peakFind_neg,~] = findpeaks(ECoG_stackN,dbs_fs,'MinPeakDistance',0.05);
                    eco_temp = squeeze(getEpochSignal(ECoG_stackN,locs-pre,locs+post));
                    
                end
                %ECoG_peak_neg{i}{j} = ECoG_peakFind_neg;
                ECoG_peak_neg{i}{j} = max(eco_temp,[],1);
            end
            
            if strcmp(ccepAnalysis,'y')
                eco_temp = squeeze(getEpochSignal(ECoG_stackP,locs-presampsCCEP,locs+postsampsCCEP));
                
                eco_epoch_ave = mean(eco_temp((tCCEP<-0.25),:),1);
                if subtractPre
                eco_temp = eco_temp - repmat(eco_epoch_ave,size(eco_temp,1),1,1);
                end
                
                ECoG_sepCCEPinternal{i}(:,:,j) = eco_temp;
            end
        end
        
    end
    
    
    if strcmp(saveOutputInternal,'y')
        
        save(fullfile(OUTPUT_DIR, ['stimInternal_',side,'_',numLeads,'DBS_' num2str(stim_chan1),'_',num2str(stim_chan2),'_fs_',num2str(stimFreq)]),...
            'ucondition','ECOG_fs','dbs_fs','DBS_sepCCEPinternal','ECoG_sepCCEPinternal','tCCEP','presamps','postsamps','pre','post');
        
    end
    
end
%% average plots of internal CCEP
if strcmp(plotCCEP,'y')
    
    prompt = {'What is the condition of interest?'};
    dlg_title = 'Condition of interest ';
    num_lines = 1;
    defaultans = {'4'};
    answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
    
    condOfInt = str2num(answer{1});
    
    % mean subtract
    p = numSubplots(numEco);

    %figure
    for j = 1:numEco
    subplot(p(1),p(2),j)
        %
        % tempEco = squeeze(ECoG_sepCCEPinternal{condOfInt}(1:length(tCCEP),:,j));
        %tempEcoBase = mean(tempEco(tCCEP<-0.25,:),1);
        %tempEcoNormalized = tempEco - repmat(tempEcoBase,[size(tempEco,1),1]);
        %mu = mean(tempEcoNormalized,2);
        %stdError = std(tempEcoNormalized,[],2)/sqrt(size(tempEcoNormalized,2));
        
        mu = mean(squeeze(ECoG_sepCCEPinternal{condOfInt}(1:length(tCCEP),:,j)),2);
        stdError = std(squeeze(ECoG_sepCCEPinternal{condOfInt}(1:length(tCCEP),:,j)),[],2)/sqrt(size(squeeze(ECoG_sepCCEPinternal{condOfInt}(1:length(tCCEP),:,j)),2));
        plot(tCCEP,mu)
        hold on
        plot(tCCEP, mu+stdError, ':');
        hold on;
        plot(tCCEP, mu-stdError, ':');
        ylabel('Voltage (V)')
        xlabel('time (ms)')
        ylim([-1e-4 1e-4])
        % xlim([min(tCCEP) 5])
        
        title(['Channel ',num2str(j)])
        
    end
    %subtitle(['ECoG CCEP responses within train for condition = ' num2str(condOfInt)])
    
    figure
    p = numSubplots(numDBS);

    for j = 1:numDBS
        
    subplot(p(1),p(2),j)        
        % tempDbs = squeeze(DBS_sepCCEPinternal{condOfInt}(1:length(tCCEP),:,j));
        % tempDbsBase = mean(tempDbs(tCCEP<-0.25,:),1);
        % tempDbsNormalized = tempDbs - repmat(tempDbsBase,[size(tempEco,1),1]);
        % mu = mean(tempDbsNormalized,2);
        % stdError = std(tempDbsNormalized,[],2)/sqrt(size(tempDbsNormalized,2));
        
        
        mu = mean(squeeze(DBS_sepCCEPinternal{condOfInt}(1:length(tCCEP),:,j)),2);
        stdError = std(squeeze(DBS_sepCCEPinternal{condOfInt}(1:length(tCCEP),:,j)),[],2)/sqrt(size(squeeze(DBS_sepCCEPinternal{condOfInt}(1:length(tCCEP),:,j)),2));
        
        plot(tCCEP,mu)
        hold on
        plot(tCCEP, mu+stdError, ':');
        hold on;
        plot(tCCEP, mu-stdError, ':');
        ylabel('Voltage (V)')
        xlabel('time (ms)')
        ylim([-1e-3 1e-3])
        %xlim([min(tCCEP) 5])
        
        % put a box around the stimulation channels of interest if need be
        if ismember(j,stimChans)
            ax = gca;
            ax.Box = 'on';
            ax.XColor = 'red';
            ax.YColor = 'red';
            ax.LineWidth = 2;
            title(['Channel ',num2str(j)],'color','red');
            
        else
            title(['Channel ',num2str(j)]);
        end
        
    end
    %subtitle(['DBS CCEP responses within train for condition = ' num2str(condOfInt)])
    
    return
    
end
%% not average
if strcmp(plotCCEP,'y')
    
    figure
    for j = 1:numEco
        subplot(4,4,j)
        plot(tCCEP,squeeze(ECoG_sepCCEPinternal{cond_int}(1:length(tCCEP),:,j)))
        xlabel('time (ms)')
        ylabel('Voltage (V)')
        title(['Channel ',num2str(j)])
        
    end
    
    figure
    for j = 1:numDBS
        subplot(2,4,j)
        plot(tCCEP,squeeze(DBS_sepCCEPinternal{cond_int}(1:length(tCCEP),:,j)))
        xlabel('time (ms)')
        ylabel('Voltage (V)')
        
        % put a box around the stimulation channels of interest if need be
        if ismember(j,stimChans)
            ax = gca;
            ax.Box = 'on';
            ax.XColor = 'red';
            ax.YColor = 'red';
            ax.LineWidth = 2;
            title(['Channel ',num2str(j)],'color','red');
            
        else
            title(['Channel ',num2str(j)]);
        end
        
    end
    
end
%% look at cceps inside of stimulation window
%
% % neeed ucondition and condition
%
% if strcmp(ccepAnalysis,'y')
%
%     ECoG_sepCCEPinternal = {};
%     DBS_sepCCEPinternal = {};
%
%     ecoEachStim = {};
%     dbsEachStim = {};
%
%     preCCEP = floor(0 * stim_fs); % pre time in sec
%
%     % 6
%     if stimFreq == 185
%         postCCEP = floor(6e-3 * stim_fs); % post time in sec, % modified DJC to look at up to 50 ms after
%     elseif stimFreq == 20
%         postCCEP = floor(51e-3 * stim_fs); % post time in sec, % modified DJC to look at up to 50 ms after
%     end
%
%
%
%     for i = 1:length(ucondition)
%         % stack it one after the other
%
%         ecoTemp = dataEpochedECOG(:,:,condition(stimTimes)==ucondition(i));
%         dbsTemp  = dataEpochedDBS(:,:,condition(stimTimes)==ucondition(i));
%
%         numTrials = size(ecoTemp,3);
%
%         % this is for 185 Hz condition
%         if stimFreq == 185
%             numPeaks = 92;
%         elseif stimFreq == 20
%             numPeaks = 9;
%
%         end
%
%         ECoG_sepCCEPinternal{i} = squeeze(reshape((permute(ecoTemp,[1 3 2])),[],1,16));
%         DBS_sepCCEPinternal{i} = squeeze(reshape((permute(dbsTemp,[1 3 2])),[],1,8));
%
%         % need to pick DBS channel with clear peaks
%         dbsFind = DBS_sepCCEPinternal{i};
%
%         if stimFreq == 185
%             [~,timesFromDBS] = findpeaks(dbsFind(:,chanExtract),'MinPeakDistance',244,'NPeaks',92);
%         elseif stimFreq == 20
%             [~,timesFromDBS] = findpeaks(dbsFind(:,chanExtract),'MinPeakDistance',2441,'NPeaks',9);
%         end
%
%         dataEpochedECoGccep = squeeze(getEpochSignal(ECoG_sepCCEPinternal{i},timesFromDBS-preCCEP,timesFromDBS+postCCEP));
%         dataEpochedDBSccep = squeeze(getEpochSignal(DBS_sepCCEPinternal{i},timesFromDBS-preCCEP,timesFromDBS+postCCEP));
%
%
%         ecoEachStim{i} = dataEpochedECoGccep;
%         dbsEachStim{i} = dataEpochedDBSccep;
%
%
%     end
%
%
%     figure
%     for i = 1:numEco
%         subplot(8,2,i)
%
%         title(['Channel '])
%     end
%     subtitle('ECoG CCEPs')
%
%     figure
%     for i = 1:numDBS
%         ECoG_neg = cellfun(@(x) x*factor,ECoG_sep,'un',0);
%
%     end
%     subtitle('DBS CCEPs')
%
% end

%% ICA analysis
doICA = 0;

if doICA
    scale_factor = 50;
    numComponentsSearch = 10;
    
    
    
    %      scale_factor = 1000;
    %      numComponentsSearch = 20;
    
    plotIt = 0;
    %stimChans = [9 17 50 58 ];
    meanSub = 1;
    %
    % [subtracted_sig_matrixS_I, subtracted_sig_cellS_I,recon_artifact_matrix,recon_artifact,t] = ...
    %     ica_artifact_remove_train(t_epoch,epochedCortEco,stimChans,eco_fs,scale_factor,numComponentsSearch,plotIt,chanInt,meanSub);
    
    % 4-23-2017 - changed orderPoly from 6 to 3
    orderPoly = 1;
    chanInt = 8;
    
    
    
    % %
    %         [processedSig,~,~,~,t] = ...
    %             ica_artifact_remove_train(t_epoch,epochedCortEco,stimChans,eco_fs,scale_factor,numComponentsSearch,plotIt,chanInt,meanSub,orderPoly);
    %
    [processedSig,~,~,~,t_process] = ...
        ica_artifact_remove_train_dbs(t,ECoG_temp,[],ECOG_fs,scale_factor,numComponentsSearch,plotIt,chanInt,meanSub,orderPoly);
    
end

%% look at CCEPs aoutside of stimulation window
if strcmp(ccepAnalysis,'y')
    
end

%% Voltage analysis - depracted - DJC - 5-23-2017
%
% if strcmp(voltageAnalysis,'y')
%
%     % get average peaks of waveform
%
%     DBS_peak_pos = {};
%     ECoG_peak_pos = {};
%
%     DBS_peak_neg = {};
%     ECoG_peak_neg = {};
%
%     factor = -1;
%     ECoG_neg = cellfun(@(x) x*factor,ECoG_sep,'un',0);
%     DBS_neg = cellfun(@(x) x*factor,DBS_sep,'un',0);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % 11/1/2016 - set min peak height to 2.5e-3 , was 5e-3 before
%
%     prompt = {'Which DBS electrode to use for extract CCEPs?'};
%     dlg_title = 'Electrode to extract (for bb908 use 8 for dbs_7,8,10,11 or 5/6 for dbs_9,12 ';
%     num_lines = 1;
%     defaultans = {'8'};
%     answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
%
%     chanExtract = str2num(answer{1});
%     MinPeakH = 2.5e-3;
%
%     dbs_stackP = squeeze(dbs_condP(:,chanExtract,:));
%     dbs_stackP = dbs_stackP(:);
%
%     if stimFreq == 185
%         [DBS_peakFind_pos,locs] = findpeaks(dbs_stackP,dbs_fs,'MinPeakDistance',0.005,'NPeaks',numTotal,'MinPeakHeight',MinPeakH);
%         findpeaks(dbs_stackP,dbs_fs,'MinPeakDistance',0.005,'NPeaks',numTotal,'MinPeakHeight',MinPeakH)
%     elseif stimFreq == 20
%         [DBS_peakFind_pos,locs] = findpeaks(dbs_stackP,dbs_fs,'MinPeakDistance',0.05,'NPeaks',numTotal,'MinPeakHeight',MinPeakH);
%         findpeaks(dbs_stackP,dbs_fs,'MinPeakDistance',0.05,'NPeaks',numTotal,'MinPeakHeight',MinPeakH)
%     end
%     locs = floor(locs *stim_fs);
%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     for i = 1:length(ucondition)
%
%         dbs_condP = DBS_sep{i};
%         dbs_condN = DBS_neg{i};
%         ECoG_condP = ECoG_sep{i};
%         ECoG_condN = ECoG_neg{i};
%
%
%         for j = 1:size(dbs_condP,2)
%             dbs_stackP = squeeze(dbs_condP(:,j,:));
%             dbs_stackP = dbs_stackP(:);
%             %[DBS_peakFind_pos,~] =
%             %findpeaks(dbs_stackP,dbs_fs,'MinPeakDistance',0.005);
%             findpeaks(dbs_stackP,dbs_fs,'MinPeakDistance',0.005)
%
%             DBS_peak_pos{i}{j} = DBS_peakFind_pos;
%
%             dbs_stackN = squeeze(dbs_condN(:,j,:));
%             dbs_stackN = dbs_stackN(:);
%             [DBS_peakFind_neg,~] = findpeaks(dbs_stackN,dbs_fs,'MinPeakDistance',0.005);
%             DBS_peak_neg{i}{j} = DBS_peakFind_neg;
%
%         end
%
%         for j = 1:size(ECoG_condP,2)
%
%             ECoG_stackP = squeeze(ECoG_condP(:,j,:));
%             ECoG_stackP = ECoG_stackP(:);
%             [ECoG_peakFind_pos,~] = findpeaks(ECoG_stackP,dbs_fs,'MinPeakDistance',0.005);
%             ECoG_peak_pos{i}{j} = ECoG_peakFind_pos;
%
%             ECoG_stackN = squeeze(ECoG_condN(:,j,:));
%             ECoG_stackN = ECoG_stackN(:);
%             [ECoG_peakFind_neg,~] = findpeaks(ECoG_stackN,dbs_fs,'MinPeakDistance',0.005);
%             ECoG_peak_neg{i}{j} = ECoG_peakFind_neg;
%         end
%
%     end
%
% end

%% Plot Hist
if strcmp(plotHist,'y')
    for i = 1:length(ucondition)
        
        
        dbsTotP(i) = figure;
        hold on
        dbsIndP(i) = figure;
        hold on
        legT = {};
        legI = {'Positive Peak','Negative Peak'};
        
        for j=1:size(dbs_condP,2)
            figure(dbsIndP(i));
            subplot(4,2,j)
            hold on
            histogram(DBS_peak_pos{i}{j})
            histogram(DBS_peak_neg{i}{j})
            
            
            title(['Electrode ' num2str(j)])
            
            figure(dbsTotP(i));
            ax1(i) = subplot(2,1,1);
            hold on
            histogram(DBS_peak_pos{i}{j})
            legT{end+1} = ['Electrode ' num2str(j)];
            xlabel('Voltage (V)')
            title('Positive Peak')
            
            ax2(i) = subplot(2,1,2);
            hold on
            histogram(DBS_peak_neg{i}{j})
            xlabel('Voltage (V)')
            title('Negative Peak')
            
        end
        
        figure(dbsIndP(i));
        legend(legI)
        xlabel('Voltage (V)')
        subtitle(['DBS stim pulses for Condition ', num2str(i)])
        
        
        figure(dbsTotP(i));
        legend(legT)
        subtitle(['DBS stim pulses for Condition ', num2str(i)])
        
        %     linkaxes([ax1(i) ax2(i)],'x')
        
        % ECoG
        
        ecogTotP(i) = figure;
        hold on
        ecogIndP(i) = figure;
        hold on
        legT = {};
        legI = {'Positive Peak','Negative Peak'};
        
        for j=1:size(ECoG_condP,2)
            figure(ecogIndP(i));
            subplot(4,4,j)
            hold on
            histogram(ECoG_peak_pos{i}{j})
            histogram(ECoG_peak_neg{i}{j})
            
            
            title(['Electrode ' num2str(j)])
            
            figure(ecogTotP(i));
            ax3(i) = subplot(2,1,1);
            hold on
            histogram(ECoG_peak_pos{i}{j})
            legT{end+1} = ['Electrode ' num2str(j)];
            xlabel('Voltage (V)')
            title('Positive Peak')
            
            ax4(i) = subplot(2,1,2);
            hold on
            histogram(ECoG_peak_neg{i}{j})
            xlabel('Voltage (V)')
            title('Negative Peak')
            
        end
        
        figure(ecogIndP(i));
        legend(legI)
        xlabel('Voltage (V)')
        subtitle(['ECoG stim pulses for Condition ', num2str(i)])
        
        
        figure(ecogTotP(i));
        legend(legT)
        subtitle(['ECoG stim pulses for Condition ', num2str(i)])
        
        %     linkaxes([ax3(i) ax4(i)],'x')
        
    end
    
end


