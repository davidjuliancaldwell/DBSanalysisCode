%% DJC - 8-29-2016 - DBS analysis script 
% This is to extract the neural data 

%% initialize output and meta dir
% clear workspace - be in the directory with all scripts necessary
close all; clear all; clc
% add path for scripts to work with data tanks
addpath('./scripts')

% set path
Z_ConstantsDBS

% subject directory, change as needed
% for David 
SUB_DIR = fullfile(myGetenv('subject_dir'));

%% load in subject

% this is from my z_constants

sid = SIDS{1};

% load in tank
if (strcmp(sid, 'bb908'))

    structureData = promptForTDTrecording;
    Sing = structureData.Sing;
    Stim = structureData.Stim;
    Valu = structureData.Valu;
    Cond = structureData.Cond;
    DBSs = structureData.DBSs;
    ECOG = structureData.ECOG;
  
    dbsElectrodes = DBSs.data;
    dbs_fs = DBSs.info.SamplingRateHz;

    ECOG = ECOG.data;
    ECOG_fs = ECOG.info.SamplingRateHz; 
    
    stimBox = Stim.data;
    stim_fs = Stim.info.SamplingRateHz;
    
    stimProgrammed = Sing.data;
    
    stimSampDeliver = Cond.data(:,1);
    condition = Cond.data(:,2);
    ttlPulse = Cond.data(:,3);
    cond_fs = Cond.info.SamplingRateHz;
   
    
end

%% plot stim
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

subtitle('Stimulation Channels')

%% Sing is wave to be delivered 

% build a burst table with the timing of stimuli from the stim file 
bursts = [];

Sing1 = stimProgrammed(:,1);


% trying something like A_BuildStimTables from BetaStim


stimMask = stimSampDeliver~=0;


% sample length of train - 500 ms 
sampsEnd = floor(0.5*stim_fs);

bursts(2,:) = find(stimMask==1);
bursts(3,:) = bursts(2,:) + repmat(sampsEnd,size(bursts(2,:)));

stims = squeeze(getEpochSignal(Sing1,(bursts(2,:)-1),(bursts(3,:))+1));
t = (0:size(stims,1)-1)/stim_fs;
t = t*1e3;
figure
plot(t,stims)
xlabel('Time (ms)');
ylabel('Voltage to be delivered')
title('Voltage to be delivered')

%delay loks to be 0.2867 ms from below.

%% Plot stims with info from above

stim1 = stimBox(:,1);
stim1Epoched = squeeze(getEpochSignal(stim1,(bursts(2,:)-1),(bursts(3,:))+1));
t = (0:size(stim1Epoched,1)-1)/stim_fs;
t = t*1e3;
figure
plot(t,stim1Epoched)
xlabel('Time (ms');
ylabel('Voltage (V)');
title('Finding the delay between current output and stim delivery')

hold on

plot(t,stims)

% get the delay in stim times

delay = round(0.1434*stim_fs/1e3);

% plot the appropriately delayed signal
figure
stimTimesBegin = bursts(2,:)-1+delay;
stimTimesEnd = bursts(3,:)-1+delay;
stim1Epoched = squeeze(getEpochSignal(stim1,stimTimesBegin,stimTimesEnd));
t = (0:size(stim1Epoched,1)-1)/stim_fs;
t = t*1e3;
figure
plot(t,stim1Epoched)
hold on
plot(t,stims(1:end-2,:))
xlabel('Time (ms');
ylabel('Voltage (V)');
title('Stim voltage monitoring with delay added in')



%% extract data - NO NEED FOR THIS SINCE DATA IS SAMPLED HIGH HERE 
% 
% try and account for delay for the stim times
stimTimes = bursts(2,:)-1+delay;

% DJC 7-7-2016, changed presamps and post samps to 1 second
presamps = round(1 * stim_fs); % pre time in sec
postsamps = round(1 * stim_fs); % post time in sec, % modified DJC to look at up to 300 ms after


%% Interpolation from miah's code - NOT NOW

% uncomment this if wanting to interpolate, broken right now

% for i = 1:size(data,2)
%     presamps = round(0.025 * fs_data); % pre time in sec
%     postsamps = round(0.125 * fs_data); % post time in sec, % modified DJC to look at up to 300 ms after
%     eco = data(:,i);
%
%     edd = zeros(size(sts));
%     efs = fs_data;
%
%     temp = squeeze(getEpochSignal(eco, sts-presamps, sts+postsamps+1));
%     foo = mean(temp,2);
%     lastsample = round(0.040 * efs);
%     foo(lastsample:end) = foo(lastsample-1);
%
%     last = find(abs(zscore(foo))>1,1,'last');
%     last2 = find(abs(diff(foo))>30e-6,1,'last')+1;
%
%     zc = false;
%
%     if (isempty(last2))
%         if (isempty(last))
%             error ('something seems wrong in the triggered average');
%         else
%             ct = last;
%         end
%     else
%         if (isempty(last))
%             ct = last2;
%         else
%             ct = max(last, last2);
%         end
%     end
%
%     while (~zc && ct <= length(foo))
%         zc = sign(foo(ct-1)) ~= sign(foo(ct));
%         ct = ct + 1;
%     end
%
%     if (ct > max(last, last2) + 0.10 * efs) % marched along more than 10 msec, probably gone to far
%         ct = max(last, last2);
%     end
%
%     % DJC - 8-31-2015 - i believe this is messing with the resizing
%     % in the figures
%     %             subplot(8,8,chan);
%     %             plot(foo);
%     %             vline(ct);
%     %
%
%     % comment this part out for no interpolation
%     for sti = 1:length(sts)
%         win = (sts(sti)-presamps):(sts(sti)+postsamps+1);
%
%         % interpolation approach
%         eco(win(presamps:(ct-1))) = interp1([presamps-1 ct], eco(win([presamps-1 ct])), presamps:(ct-1));
%     end
%
%     data(:,i) = eco;
% end

%% get the data epochs


dataEpochedECOG = squeeze(getEpochSignal(ECOG,stimTimes-presamps,stimTimes+postsamps+1));
dataEpochedDBS = squeeze(getEpochSignal(dbsElectrodes,stimTimes-presamps,stimTimes+postsamps+1));

% set the time vector to be set by the pre and post samps
t = (-presamps:postsamps)*1e3/fs_data;



%% 6-23-2016
% plot aggregate data for each stim type - THIS ONLY WORKS IF THERE ARE 30
% TOTAL STIM EPOCHS - stim 9 for instance only has 28 total

% chunk out data

% to separate out low and high

%k=1:8;
%j = 1:8;

k = 1:10;
j = 1:10;

% figure
for i = 1:64
%     hold on
%     subplot(8,8,i)
%     plot(t,1e6*mean(dataEpoched(:,i,j),3),'m','LineWidth',1)
%     
%     ylim([-150 150])
%     xlim([-100 300])
%     title(sprintf('Channel %d',i))
%     %     pause(1)
    
    dataEpochedLow(:,i,k) = dataEpoched(:,i,j);
    
end
% subtitle('Average traces for all stimulations - means not subtracted - stims 1:10')
% xlabel('Time (ms)')
% ylabel('Voltage (uV)')
% 
% figure
for i = 65:80
%     hold on
%     subplot(8,8,i-64)
%     plot(t,1e6*mean(dataEpoched(:,i,j),3),'m','LineWidth',1)
%     
%     ylim([-150 150])
%     xlim([-100 300])
%     title(sprintf('Channel %d',i))
%     %     pause(1)
    dataEpochedLow(:,i,k) = dataEpoched(:,i,j);
    
    
end
% subtitle('Average traces for all stimulations - means not subtracted - stims 1:10')
% xlabel('Time (ms)')
% ylabel('Voltage (uV)')

k= 1:10;
%j = 9:18;
j = 11:20;

% figure
for i = 1:64
%     hold on
%     subplot(8,8,i)
%     plot(t,1e6*mean(dataEpoched(:,i,j),3),'m','LineWidth',1)
%     
%     ylim([-150 150])
%     xlim([-100 300])
%     title(sprintf('Channel %d',i))
%     %     pause(1)
    dataEpochedMid(:,i,k) = dataEpoched(:,i,j);
    
    
end
% subtitle('Average traces for all stimulations - means not subtracted stims - stims 11:20')
% xlabel('Time (ms)')
% ylabel('Amplitude (\muV)')
% 
% figure
for i = 65:80
%     hold on
%     subplot(8,8,i-64)
%     plot(t,1e6*mean(dataEpoched(:,i,j),3),'m','LineWidth',1)
%     
%     ylim([-150 150])
%     xlim([-100 300])
%     title(sprintf('Channel %d',i))
%     %     pause(1)
    dataEpochedMid(:,i,k) = dataEpoched(:,i,j);
    
    
end
% subtitle('Average traces for all stimulations - means not subtracted - stims 11:20 ')
% xlabel('Time (ms)')
% ylabel('Amplitude (\muV)')

%k=1:10;
j = 21:30;
%j = 21:39;

% figure
for i = 1:64
%     hold on
%     subplot(8,8,i)
%     plot(t,1e6*mean(dataEpoched(:,i,j),3),'m','LineWidth',1)
%     
%     ylim([-150 150])
%     xlim([-100 300])
%     title(sprintf('Channel %d',i))
%     %     pause(1)
    dataEpochedHigh(:,i,k) = dataEpoched(:,i,j);
    
    
end
% subtitle('Average traces for all stimulations - means not subtracted - stims 21:30 ' )
% xlabel('Time (ms)')
% ylabel('Amplitude (\muV)')
% 
% figure
for i = 65:80
%     hold on
%     subplot(8,8,i-64)
%     plot(t,1e6*mean(dataEpoched(:,i,j),3),'m','LineWidth',1)
%     
%     ylim([-150 150])
%     xlim([-100 300])
%     title(sprintf('Channel %d',i))
%     %     pause(1)
    dataEpochedHigh(:,i,k) = dataEpoched(:,i,j);
    
    
end
% subtitle('Average traces for all stimulations - means not subtracted - stims 21:30')
% xlabel('Time (ms)')
% ylabel('Amplitude (\muV)')


%% 6-23-2016 - plot channel of interest
% 
% % pick channel
% i = 21;
% % pick range of stims
% j = 1:10;
% %j = 11:20;
% %j = 21:30;
% 
% figure
% plot(t,1e6*mean(dataEpoched(:,i,j),3),'m','LineWidth',1)
% xlabel('time (ms)')
% ylabel('Amplitude (\muV)')
% title(['Average for subselected stims for channel ', num2str(i)])

