%% 10-24-2016 - Quick DBS stim extraction - David Caldwell - script to look at stim spacing
% requires getEpochSignal.m , subtitle.m , numSubplots.m , vline.m
% djc - 2/8/2018 for 5e0cf
% 7.12.2018
% David.J.Caldwell
%% initialize output and meta dir
% clear workspace
close all; clear all; clc
SIDS = {'5e0cf','b26b7','80301','46c2'};
OUTPUT_DIR = pwd;

%%

% load in the datafile of interest!
% have to have a value assigned to the file to have it wait to finish
% loading...mathworks bug
sid = SIDS{4};
[structureData,filepath] = promptForTDTrecording;
Sing = structureData.Sing;
Stim = structureData.Stim;
DBSs = structureData.DBSs;
ECOG = structureData.ECOG;

filenameCell = strsplit(filepath,'\');
fileName_mat = strsplit(filenameCell{end},'.');
fileName = fileName_mat{1};

dbsElectrodes = DBSs.data;
ECOGelectrodes = ECOG.data;

switch sid
    case '5e0cf'
    dbsElectrodes = dbsElectrodes(:,1:4);
        ECOGelectrodes = ECOGelectrodes(:,1:8);
    case 'b26b7'
           dbsElectrodes = dbsElectrodes(:,1:8);
        ECOGelectrodes = ECOGelectrodes(:,1:8); 
    case '46c2a'
                   dbsElectrodes = [];
        ECOGelectrodes = ECOGelectrodes(:,1:8); 
end

ECOG_fs = ECOG.info.SamplingRateHz;
dbs_fs = DBSs.info.SamplingRateHz;

stimBox = Stim.data;
stim_fs = Stim.info.SamplingRateHz;

stimProgrammed = Sing.data;

%%
% ui box for input for stimulation channels
prompt = {'how many channels did we record from? e.g 8 ', 'what were the stimulation channels? e.g 7 8 ', 'how long before each stimulation do you want to look? in ms e.g. 1', 'how long after each stimulation do you want to look? in ms e.g 5','process data to remove z>3 outliers?','subtract mean if DC coupled?','bad channels?','save plots?'};
dlg_title = 'StimChans';
num_lines = 1;
defaultans = {'16','5 6','1','5','n','n','','y'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
numChans = str2num(answer{1});
chans = str2num(answer{2});
preTime = str2num(answer{3});
postTime = str2num(answer{4});
zScoreThresh = answer{5};
meanSubtract = answer{6};
savePlot = answer{8};
badChans = str2num(answer{7});

%%
% first and second stimulation channel
stim_1 = chans(1);
stim_2 = chans(2);
stimChans = [stim_1 stim_2];
elec_vec = [1:numChans];
badChans = [stimChans badChans];

% get sampling rates
fs_data = DBSs.info.SamplingRateHz;
fs_stim = Stim.info.SamplingRateHz;

% stim data
stim = Stim.data;

% current data
sing = Sing.data;

% recording data - 5-23-2017, DJC, add in ECoG
%data = DBSs.data;

data = [ECOGelectrodes,dbsElectrodes, ];

data = data(:,1:8);

% DJC - 8-24-2016 - if it's DC coupled, below

% zscore threshold to clear it

if strcmp(zScoreThresh,'y')
    zScoredSig = zscore(data(:,1));
    dataTemp = data(abs(zScoredSig)<3,:);
    clear data;
    data = dataTemp;
    clear dataTemp;
end

if strcmp(meanSubtract,'y')
    data = data-repmat(mean(data,1), [size(data, 1), 1]);
end

%% plot stim channels if interested

prompt = {'plot intermediate figures that show stim voltage, currents, and delays? "y" or "n" '};
dlg_title = 'StimChans';
num_lines = 1;
defaultans = {'y'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
plotIt = answer(1);

if strcmp(plotIt,'y')
    figure;
    hold on;
    for i = 1:size(stim,2)
        t = (0:length(stim)-1)/fs_stim;
        subplot(2,2,i);
        plot(t*1e3,stim(:,i));
        title(sprintf('Channel %d',i))
    end
    
    xlabel('Time (ms)');
    ylabel('Amplitude (V)');
    
    %subtitle('Stimulation Channels');
end

%% Sing looks like the wave to be delivered, with amplitude in uA


% build a burst table with the timing of stimuli
bursts = [];

% first channel of current
Sing1 = sing(:,1);
fs_sing = Sing.info.SamplingRateHz;

samplesOfPulse = round(2*fs_stim/1e3);

Sing1Mask = Sing1~=0;
dmode = diff([0 Sing1Mask' 0 ]);

dmode(end-1) = dmode(end);

bursts(2,:) = find(dmode==1);
bursts(3,:) = find(dmode==-1);

singEpoched = squeeze(getEpochSignal(Sing1,(bursts(2,:)-1),(bursts(3,:))+1));
t = (0:size(singEpoched,1)-1)/fs_sing;
t = t*1e3;

if strcmp(plotIt,'y')
    
    figure
    plot(t,singEpoched)
    xlabel('Time (ms)');
    ylabel('Current to be delivered (\muA)')
    title('Current to be delivered for all trials')
end


%% Plot stims with info from above, and find the delay!

stim1stChan = stim(:,1);
stim1Epoched = squeeze(getEpochSignal(stim1stChan,(bursts(2,:)-1),(bursts(3,:))+120));

% put around 0
stim1Epoched = stim1Epoched - repmat(stim1Epoched(1,:),size(stim1Epoched,1),1);
t = (0:size(stim1Epoched,1)-1)/fs_stim;
t = t*1e3;

% delay looks to be 7 samples

%% DJC - 10-28-2016 - normalize to baseline

labels = max(singEpoched);

uniqueLabels = unique(labels);
%
% figure
%
% if strcmp(plotIt,'y')
%     for i = uniqueLabels
%         stimInterest = stim1Epoched(:,labels==i);
%         baseline = mean(stim1Epoched(1:5,:));
%         baselineRepped = repmat(baseline,[size(stim1Epoched,1), 1]);
%         stimNorm = stimInterest-baselineRepped(:,labels==i);
%         plot(t,stimNorm);
%         hold on
%     end
%     xlabel('Time (ms)');
%     ylabel('Voltage (V)');
%     title('Individual Stim Currents overlaid')
%
% end
%%
% legLabels = {[num2str(uniqueLabels(1))]};
%
% k = 2;
% if length(uniqueLabels>1)
%     for i = uniqueLabels(2:end)
%         legLabels{end+1} = [num2str(uniqueLabels(k))];
%         k = k+1;
%     end
% end
%
% legend(legLabels);

% plot divide by current
%
% figure
% stimDivideTotal = [];
% if strcmp(plotIt,'y')
%     for i = uniqueLabels(uniqueLabels~=1)
%         stimInterest = stim1Epoched(:,labels==i);
%         baseline = mean(stim1Epoched(1:5,:));
%         baselineRepped = repmat(baseline,[size(stim1Epoched,1), 1]);
%         stimNorm = stimInterest-baselineRepped(:,labels==i);
%         stimDivide = stimNorm./(i*1e-6);
%         plot(t,stimDivide);
%         stimDivideTotal = [stimDivideTotal stimDivide];
%         hold on
%     end
%     xlabel('Time (ms)');
%     ylabel('V/I');
%     title('Voltage divided by current')
%
% end

%%


% if strcmp(plotIt,'y')
%     figure
%     plot(t,stim1Epoched)
%     xlabel('Time (ms)');
%     ylabel('Voltage (V)');
%     title('Finding the delay between current output and stim delivery')
% end


% get the delay in stim times - looks to be 7 samples or so
delay = round(0.1434*fs_stim/1e3);

delay = 0; %%%% setting delay = 0 to show better plots

% plot the appropriately delayed signal
if strcmp(plotIt,'y')
    stimTimesBegin = bursts(2,:)-1+delay;
    stimTimesEnd = bursts(3,:)-1+delay+120;
    % stim1Epoched = squeeze(getEpochSignal(stim1stChan,stimTimesBegin,stimTimesEnd));
    t = (0:size(stim1Epoched,1)-1)/fs_stim;
    t = t*1e3;
    figure
    plot(t,stim1Epoched(:,:))
    xlabel('Time (ms)');
    ylabel('Voltage (V)');
    title('Stim voltage')
    
    if savePlot
        SaveFig(OUTPUT_DIR,['stimChans_' num2str(stimChans(1)) '_' num2str(stimChans(2)) '_' 'voltageMonitor']);
    end
end


% % redefine delay for other work - 10-28-2016
%
% delay = round(0.1434*fs_stim/1e3);


%% extract data

% try and account for delay for the stim times
stimTimes = bursts(2,:)-1+delay;

% DJC 7-7-2016, changed presamps and postsamps to be user defined
presamps = round(preTime/1000 * fs_data); % pre time in sec
postsamps = round(postTime/1000 * fs_data); % post time in sec,


% sampling rate conversion between stim and data
fac = fs_stim/fs_data;

% find times where stims start in terms of data sampling rate
sts = round(stimTimes / fac);


% looks like there's an additional 14 sample delay between the stimulation being set to
% be delivered....and the ECoG recording. which would be 2.3 ms?

delay2 = 14;
sts = round(stimTimes / fac) + delay2;
%sts = round(stimTimes / fac);

%% get the data epochs

%%%%%%%%%%%%%%%%%%%%%% for stim param 12 at first pass

%data = data(1:6e6,:);
%sts = sts(sts<6e6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dataEpoched = squeeze(getEpochSignal(data,sts-presamps,sts+postsamps+1));

% set the time vector to be set by the pre and post samps
t = (-presamps:postsamps)*1e3/fs_data;

% get rid of bad channels
chansVec_goods = ones(numChans,1);
chansVec_goods(badChans) = 0;
dataEpoched(:,~chansVec_goods,:) = nan;


%% make the decision to scale it

% ui box for input
prompt = {'scale the y axis to the maximum stim pulse value? "y" or "n" '};
dlg_title = 'Scale';
num_lines = 1;
defaultans = {'n'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
scaling = answer{1};


if strcmp(scaling,'y')
    maxVal = max(dataEpoched(:));
    minVal = min(dataEpoched(:));
end


%% plot individual trials for each condition on a different graph

labels = max(singEpoched);

%%%%%%%%%%%%%%%%%%%% for stim param 12 at first pass
%labels = labels(sts<6e6);
%%%%%%%%%%%%%%%%%%%%

uniqueLabels = unique(labels);

% intialize counter for plotting
k = 1;

% make vector of stim channels
% determine number of subplot
subPlots = numSubplots(numChans);
p = subPlots(1);
q = subPlots(2);

% plot each condition separately e.g. 1000 uA, 2000 uA, and so on

for i=uniqueLabels
    figure('units','normalized','outerposition',[0 0 1 1]);
    dataInterest = dataEpoched(:,:,labels==i);
    for j = 1:numChans
        subplot(p,q,j);
        plot(t,squeeze(dataInterest(:,j,:)));
        xlim([min(t) max(t)]);
        
        % change y axis scaling if necessary
        if strcmp(scaling,'y')
            ylim([minVal maxVal]);
        end
        
        % put a box around the stimulation channels of interest if need be
        if ismember(j,stimChans)
            ax = gca;
            ax.Box = 'on';
            ax.XColor = 'red';
            ax.YColor = 'red';
            ax.LineWidth = 2;
            title(num2str(j),'color','red');
            
        else
            title(num2str(j));
            
        end
        vline(0);
        
    end
    
    % label axis
    xlabel('time in ms');
    ylabel('voltage in \muV');
    %subtitle(['Individual traces - Current set to ',num2str(uniqueLabels(k)),' \muA']);
    % subtitle(['Individual traces - Voltage set to ',num2str(uniqueLabels(k)),' V']);
    
    % get cell of raw values, can use this to analyze later
    dataEpochedCell{k} = dataInterest;
    
    % get averages to plot against each for later
    % cell function, can use this to analyze later
    dataAvgs{k} = mean(dataInterest,3);
    %increment counter
    k = k + 1;
    
    if savePlot
        SaveFig(OUTPUT_DIR,['stimChans_' num2str(stimChans(1)) '_' num2str(stimChans(2)) '_' 'individualRecordings']);
    end
    
end

%% plot averages for 3 conditions on the same graph

k = 1;
figure('units','normalized','outerposition',[0 0 1 1]);
for k = 1:length(dataAvgs)
    
    tempData = dataAvgs{k};
    
    for j = 1:numChans
        s = subplot(p,q,j);
        plot(t,squeeze(tempData(:,j)),'linewidth',2);
        hold on;
        xlim([min(t) max(t)]);
        
        
        % change y axis scaling if necessary
        if strcmp(scaling,'y')
            ylim([minVal maxVal]);
        end
        
        
        if ismember(j,stimChans)
            ax = gca;
            ax.Box = 'on';
            ax.XColor = 'red';
            ax.YColor = 'red';
            ax.LineWidth = 2;
            title(num2str(j),'color','red')
            
        else
            title(num2str(j));
            
        end
        
        vline(0);
        
    end
    gcf;
end
xlabel('time in ms');
ylabel('voltage in \muV');
%subtitle(['Averages for all conditions']);
legLabels = {[num2str(uniqueLabels(1))]};

k = 2;
if length(uniqueLabels>1)
    for i = uniqueLabels(2:end)
        legLabels{end+1} = [num2str(uniqueLabels(k))];
        k = k+1;
    end
    
end

legend(s,legLabels);


if savePlot
    SaveFig(OUTPUT_DIR,['stimChans_' num2str(stimChans(1)) '_' num2str(stimChans(2)) '_' 'averageRecordings']);
end

%% histograms
elec_vec = [1:size(dataEpoched,2)];
middlePts_1st = squeeze(dataEpoched(83,:,:));
middlePts_2nd = squeeze(dataEpoched(139,:,:));

mean_1st = abs(mean(middlePts_1st,2));
mean_2nd = abs(mean(middlePts_2nd,2));
std_1st = std(middlePts_1st,[],2);
std_2nd = std(middlePts_2nd,[],2);

figure
hold on
errorbar(elec_vec,mean_1st,std_1st,'o')
errorbar(elec_vec,mean_2nd,std_2nd,'o')
xlabel('Electrode')
ylabel('\muV')
title(['Mean and std for middle of recorded pulses - stim chans ' num2str(stimChans(1)) ' _ ' num2str(stimChans(2))])
%vline(9,'k')
vline(stimChans(1),'g')
vline(stimChans(2),'b')
legend('1st phase','2nd phase','DBS/ECoG','- chan','+ chan')
%%
if savePlot
    SaveFig(OUTPUT_DIR,['stimChans_' num2str(stimChans(1)) '_' num2str(stimChans(2)) '_' 'meanStd']);
end

%% save it - djc 2/8/2018
saveData = 1;
if saveData
    OUTPUT_DIR = pwd;
    fs = fs_data;
    save(sprintf(['stimSpacingDBS-%s-stim_%d-%d'], sid, stimChans(1),stimChans(2)),'dataEpoched','dataEpochedCell','stimChans','t','singEpoched','stim1Epoched','middlePts_1st','middlePts_2nd','fs');
end