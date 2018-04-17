%% 3-16-2018 - script to extract ECoG/DBS recordings during Cadwell recordings
% of note is that we don't have clear indication of when cortical stim was
% delivered, so will have to use threshold crossing or another technique
%close all;clear all;clc
DATA_DIR = 'G:\My Drive\GRIDLabDavidShared\DBS';
% set path, set
Z_ConstantsDBS_cadwell

% subject directory, change as needed
% for David
SUB_DIR = fullfile(myGetenv('subject_dir'));

%% load in subject

% this is from my z_constants
%SIDS = {'80301','63ce7','1dd75','56a68','b305e','329c6','5e0cf','c1c8c','b26b7'};

% for param sweep, look at subjects 1,2,9
%sid = input('what is the sid?\n','s');
%sid = SIDS{2}; % MUST SWITCH THIS, either 1,2,9
sid = '50ad9';
% load in tank
switch sid
        
    case '80301'
        [structureData,filepath] = promptForTDTrecording(DATA_DIR);
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
        
    case '63ce7'
        [structureData,filepath] = promptForTDTrecording(DATA_DIR);
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
    case '1dd75'
        [structureData,filepath] = promptForTDTrecording(DATA_DIR);
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
        ttlPulse = Cond.data(:,3);
        cond_fs = Cond.info.SamplingRateHz;
        
    otherwise
        [structureData,filepath] = promptForTDTrecording(DATA_DIR);
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
%dbsStim = [2 3] + 1;
dbsStim = [];
ECoGStim = [5];

% plot stim
%
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


%
ECoGdata = ECOG.data;
DBSdata = DBSs.data;
    t = 1e3*(0:length(stimBox)-1)/ECOG_fs;


ECoGdata(:,ECoGStim) = 0;
DBSdata(:,dbsStim) = 0;

figure
numEco = size(ECoGdata,2);
p = numSubplots(numEco);
for j = 1:numEco
    axEeco(j) = subplot(p(1),p(2),j);
    plot(t,ECoGdata(:,j));
    title(['ECoG Channel ' num2str(j)]);
end
xlabel('time (ms)')
ylabel('voltage (V)')
%subtitle(['ECoG Electrodes, Condition ' num2str(i)]);
linkaxes(axEeco,'xy');


% plot DBS Electrodes

figure
numDbs = size(DBSdata,2);
p = numSubplots(numDbs);
for j = 1:numDbs
   axDbs(j)= subplot(p(1),p(2),j);
    plot(t,DBSdata(:,j));
    title(['DBS Channel ' num2str(j)]);
end
xlabel('time (ms)')
ylabel('voltage (V)')
%  subtitle(['DBS Electrodes, Condition ' num2str(i)]);
linkaxes(axDbs,'xy');
%
ECoGchanInt = 6;
DBSchanInt = 2;

figure
plot(t,ECoGdata(:,ECoGchanInt));
title(['ECoG channel ' num2str(ECoGchanInt)])
xlabel('time (ms)')
ylabel('voltage (V)')

figure
plot(t,DBSdata(:,DBSchanInt));
title(['DBS channel ' num2str(DBSchanInt)])
xlabel('time (ms)')
ylabel('voltage (V)')


