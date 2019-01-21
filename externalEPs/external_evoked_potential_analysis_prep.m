%% DJC - 3-29-2018
% This is to extract the neural data outside of the DBS trains to look at
% evoked potentials

%% initialize output and meta dir
% clear workspace - be in the directory with all scripts necessary
close all; clear all; clc

% set path, set
Z_ConstantsDBS_externalEPs

% subject directory, change as needed
% for David
SUB_DIR = fullfile(myGetenv('subject_dir'));

%% load in subject

% here are subjects that we have acquired MEP data on
%SIDS = {'80301','63ce7','1dd75','56a68','b305e','329c6','c1c8c','b26b7'};


sid = 'c1c8c';
% load in tank
switch sid
    case 'bb908' % Gpi pateitn
        
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
    case '63ce7'
        % no paramsweep
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
        % we have to reverse the order of the DBS electrodes
        dbsElectrodes = fliplr(dbsElectrodes(:,1:4));
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
    case '56a68'
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
    case 'c1c8c'
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
    case '5e0cf'
        % no paramsweep
    case '329c6'
        % no paramsweep
    case 'b305e'
                % paramsweep one side
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
    case 'b26b7'
        % paramsweep one side
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
    'Plot time series of DBS and ECoG electrodes? y or n',
    'Plot Specific channels or conditions of interest? y or n',...
    'Save output file'};
dlg_title = 'Channel of Interest';
num_lines = 1;
defaultans = {'y','y','y','y','y','y'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

plotStim = answer{1};
plotTime = answer{2};
plotCond = answer{3};
saveOutput = answer{4};

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
defaultans = {'0.1','0.795','180'};
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
defaultans = {'L','single','1','2'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

side = answer{1};
numLeads = answer{2};
stim_chan1 = str2num(answer{3});
stim_chan2 = str2num(answer{4});
stimChans = [stim_chan1 stim_chan2];

if strcmp(saveOutput,'y')
    
    
    save(fullfile(OUTPUT_DIR, [sid '_stim_',side,'_',numLeads,'DBS_' num2str(stim_chan1),'_',num2str(stim_chan2)]),...
        'ucondition','ECOG_fs','dbs_fs','DBS_sep','ECoG_sep','t','presamps','postsamps','pre','post','stimChans');
    
end

windowInt = [-500 2500];
plot_EPs(ucondition,DBS_sep,ECoG_sep,t,stimChans,windowInt)