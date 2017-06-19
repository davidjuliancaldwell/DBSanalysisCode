%% DJC 11-2-2016 - updated 6-6-2017
% Script to subtract opposite stimulation waveform pairs
close all;clear all;clc
Z_ConstantsDBS


load(fullfile(OUTPUT_DIR,'stimInternal_l_singleDBS_3_4_fs_185.mat'));
%load(fullfile(OUTPUT_DIR,'stimInternal_R_bothDBS_7_8_fs_185.mat'));

ECoG_sepCCEPinternal1 = ECoG_sepCCEPinternal;
DBS_sepCCEPinternal1 = DBS_sepCCEPinternal;


load(fullfile(OUTPUT_DIR,'stimInternal_l_singleDBS_4_3_fs_185.mat'));
%load(fullfile(OUTPUT_DIR,'stimInternal_R_bothDBS_8_7_fs_185.mat'));

ECoG_sepCCEPinternal2 = ECoG_sepCCEPinternal;
DBS_sepCCEPinternal2 = DBS_sepCCEPinternal;
%%

prompt = {'What is the condition of interest e.g. 4 ?','what are the stim channels? e.g. [1 2]',...
    'what are the channels to exclude from rerereferencing e.g [1 2 3]','"median", "mean","singleChan","bipolar", or "n" rereference for ECoG','"median", "mean", "singleChan","bipolar", or "n" rereference for DBS'};
dlg_title = 'Condition of interest ';
num_lines = 1;
defaultans = {'4','[1 2]','[9 10 11 12]','mean','n'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

condOfInt = str2num(answer{1});
stimChans = str2num(answer{2});
badChans = str2num(answer{3});
rrEco = answer{4};
rrDbs = answer{5};

numEco = size(ECoG_sepCCEPinternal{1},3);
numDBS = size(DBS_sepCCEPinternal{1},3);

%% rereference
switch rrEco
    case 'median'
        for i = 1:length(ECoG_sepCCEPinternal)
            
            %%%%%%%%%%%%%%%%%%%%%%%
            tempEco1 = squeeze(ECoG_sepCCEPinternal1{condOfInt}(1:length(tCCEP),:,:));
            tempEcoNormalized = rereference_CAR_median(tempEco1,rrEco,badChans,[1 3 2]);
            ECoG_sepCCEPinternal1{condOfInt}(1:length(tCCEP),:,:) = tempEcoNormalized;
            %%%%%%%%%%%%%%%%%%%%%%%
            tempEco2 = squeeze(ECoG_sepCCEPinternal2{condOfInt}(1:length(tCCEP),:,:));
            tempEcoNormalized = rereference_CAR_median(tempEco2,rrEco,badChans,[1 3 2]);
            ECoG_sepCCEPinternal2{condOfInt}(1:length(tCCEP),:,:) = tempEcoNormalized;
            
        end
    case 'mean'
        
        %%%%%%%%%%%%%%%%%%%%%%%
        tempEco1 = squeeze(ECoG_sepCCEPinternal1{condOfInt}(1:length(tCCEP),:,:));
        tempEcoNormalized = rereference_CAR_median(tempEco1,rrEco,badChans,[1 3 2]);
        ECoG_sepCCEPinternal1{condOfInt}(1:length(tCCEP),:,:) = tempEcoNormalized;
        %%%%%%%%%%%%%%%%%%%%%%%
        tempEco2 = squeeze(ECoG_sepCCEPinternal2{condOfInt}(1:length(tCCEP),:,:));
        tempEcoNormalized = rereference_CAR_median(tempEco2,rrEco,badChans,[1 3 2]);
        ECoG_sepCCEPinternal2{condOfInt}(1:length(tCCEP),:,:) = tempEcoNormalized;
    case 'singleChan'
                %%%%%%%%%%%%%%%%%%%%%%%
        tempEco1 = squeeze(ECoG_sepCCEPinternal1{condOfInt}(1:length(tCCEP),:,:));
        tempEcoNormalized = rereference_CAR_median(tempEco1,rrEco,badChans,[1 3 2],1);
        ECoG_sepCCEPinternal1{condOfInt}(1:length(tCCEP),:,:) = tempEcoNormalized;
        %%%%%%%%%%%%%%%%%%%%%%%
        tempEco2 = squeeze(ECoG_sepCCEPinternal2{condOfInt}(1:length(tCCEP),:,:));
        tempEcoNormalized = rereference_CAR_median(tempEco2,rrEco,badChans,[1 3 2],1);
        ECoG_sepCCEPinternal2{condOfInt}(1:length(tCCEP),:,:) = tempEcoNormalized;
    case 'bipolar'
                %%%%%%%%%%%%%%%%%%%%%%%
        tempEco1 = squeeze(ECoG_sepCCEPinternal1{condOfInt}(1:length(tCCEP),:,:));
        tempEcoNormalized = rereference_CAR_median(tempEco1,rrEco,badChans,[1 3 2]);
        ECoG_sepCCEPinternal1{condOfInt}(1:length(tCCEP),:,:) = tempEcoNormalized;
        %%%%%%%%%%%%%%%%%%%%%%%
        tempEco2 = squeeze(ECoG_sepCCEPinternal2{condOfInt}(1:length(tCCEP),:,:));
        tempEcoNormalized = rereference_CAR_median(tempEco2,rrEco,badChans,[1 3 2]);
        ECoG_sepCCEPinternal2{condOfInt}(1:length(tCCEP),:,:) = tempEcoNormalized;
        
end

switch rrDbs
    case 'median'
        for i = 1:length(ECoG_sepCCEPinternal)
            
            
            %%%%%%%%%%%%%%%%%%%%%%%
            tempDbs1 = squeeze(DBS_sepCCEPinternal1{condOfInt}(1:length(tCCEP),:,:));
            tempDbsNormalized = rereference_CAR_median(tempDbs1,rrDbs,stimChans,[1 3 2]);
            DBS_sepCCEPinternal1{condOfInt}(1:length(tCCEP),:,:) = tempDbsNormalized;
            %%%%%%%%%%%%%%%%%%%%%%%
            tempDbs2 = squeeze(DBS_sepCCEPinternal2{condOfInt}(1:length(tCCEP),:,:));
            tempDbsNormalized = rereference_CAR_median(tempDbs2,rrDbs,stimChans,[1 3 2]);
            DBS_sepCCEPinternal2{condOfInt}(1:length(tCCEP),:,:) = tempDbsNormalized;
            
        end
    case 'mean'
        
        %%%%%%%%%%%%%%%%%%%%%%%
        tempDbs1 = squeeze(DBS_sepCCEPinternal1{condOfInt}(1:length(tCCEP),:,:));
        tempDbsNormalized = rereference_CAR_median(tempDbs1,rrDbs,stimChans,[1 3 2]);
        DBS_sepCCEPinternal1{condOfInt}(1:length(tCCEP),:,:) = tempDbsNormalized;
        %%%%%%%%%%%%%%%%%%%%%%%
        tempDbs2 = squeeze(DBS_sepCCEPinternal2{condOfInt}(1:length(tCCEP),:,:));
        tempDbsNormalized = rereference_CAR_median(tempDbs2,rrDbs,stimChans,[1 3 2]);
        DBS_sepCCEPinternal2{condOfInt}(1:length(tCCEP),:,:) = tempDbsNormalized;
    case 'singleChan'
                %%%%%%%%%%%%%%%%%%%%%%%
        tempDbs1 = squeeze(DBS_sepCCEPinternal1{condOfInt}(1:length(tCCEP),:,:));
        tempDbsNormalized = rereference_CAR_median(tempDbs1,rrDbs,stimChans,[1 3 2],1);
        DBS_sepCCEPinternal1{condOfInt}(1:length(tCCEP),:,:) = tempDbsNormalized;
        %%%%%%%%%%%%%%%%%%%%%%%
        tempDbs2 = squeeze(DBS_sepCCEPinternal2{condOfInt}(1:length(tCCEP),:,:));
        tempDbsNormalized = rereference_CAR_median(tempDbs2,rrDbs,stimChans,[1 3 2],1);
        DBS_sepCCEPinternal2{condOfInt}(1:length(tCCEP),:,:) = tempDbsNormalized;
    case 'bipolar'
                    %%%%%%%%%%%%%%%%%%%%%%%
            tempDbs1 = squeeze(DBS_sepCCEPinternal1{condOfInt}(1:length(tCCEP),:,:));
            tempDbsNormalized = rereference_CAR_median(tempDbs1,rrDbs,stimChans,[1 3 2]);
            DBS_sepCCEPinternal1{condOfInt}(1:length(tCCEP),:,:) = tempDbsNormalized;
            %%%%%%%%%%%%%%%%%%%%%%%
            tempDbs2 = squeeze(DBS_sepCCEPinternal2{condOfInt}(1:length(tCCEP),:,:));
            tempDbsNormalized = rereference_CAR_median(tempDbs2,rrDbs,stimChans,[1 3 2]);
            DBS_sepCCEPinternal2{condOfInt}(1:length(tCCEP),:,:) = tempDbsNormalized;
end
%%
% line vs. shaded

plotTypeVec = {'doShade','error','individual'};
plotType = plotTypeVec{1};

subtractPre = 0;

figure
p = numSubplots(numEco);
p(1) = 2
p(2) = 8
for j = 1:numEco
    subplot(p(1),p(2),j)
    % positive polarity
    tempEco1 = squeeze(ECoG_sepCCEPinternal1{condOfInt}(1:length(tCCEP),:,j));
    tempEcoBase1 = mean(tempEco1(tCCEP<-0.25,:),1);
    
    % negative polarity
    tempEco2 = squeeze(ECoG_sepCCEPinternal2{condOfInt}(1:length(tCCEP),:,j));
    tempEcoBase2 = mean(tempEco2(tCCEP<-0.25,:),1);
    if subtractPre
        tempEcoNormalized2 = tempEco2 - repmat(tempEcoBase2,[size(tempEco2,1),1]);
        tempEcoNormalized1 = tempEco1 - repmat(tempEcoBase1,[size(tempEco1,1),1]);
        
    else
        tempEcoNormalized2 = tempEco2;
        tempEcoNormalized1 = tempEco1;
        
    end
    
    
    % normalize if one is larger than the other
    
    if size(tempEcoNormalized2,2)>size(tempEcoNormalized1,2)
        tempEcoNormalized2 = tempEcoNormalized2(:,1:size(tempEcoNormalized1,2));
    else
        tempEcoNormalized1 = tempEcoNormalized1(:,1:size(tempEcoNormalized2,2));
    end
    
    
    tempEcoNormalized = [tempEcoNormalized1 tempEcoNormalized2];
    
    % rereference if wanted
    
    mu = mean(tempEcoNormalized,2);
    stdError = std(tempEcoNormalized,[],2)/sqrt(size(tempEcoNormalized,2));
    
    %mu = mean(squeeze(ECoG_sepCCEPinternal{condOfInt}(1:length(tCCEP),:,j)),2);
    % stdError = std(squeeze(ECoG_sepCCEPinternal{condOfInt}(1:length(tCCEP),:,j)),[],2)/sqrt(size(squeeze(ECoG_sepCCEPinternal{condOfInt}(1:length(tCCEP),:,j)),2));
    switch plotType
        case 'errorPlot'
            plot(tCCEP,mu)
            hold on
            plot(tCCEP, mu+stdError, ':');
            hold on;
            plot(tCCEP, mu-stdError, ':');
        case 'doShade'
            plotBTLError(tCCEP,tempEcoNormalized,'CI');
        case 'individual'
            plot(tCCEP,tempEcoNormalized)
            hold on
            plot(mu,'k','linewidth',2)
    end
    
    ylim([-5e-6 5e-6])
    xlim([0 5])
    
    %xlim([min(tCCEP) 5])
    
    title(['ECoG Channel ',num2str(j)])
    
end
ylabel('Voltage (V)')
xlabel('time (ms)')
figure
p = numSubplots(numDBS);
p(1) = 2
p(2) = 4
for j = 1:numDBS
    
    subplot(p(1),p(2),j)
    % one polarity
    tempDbs1 = squeeze(DBS_sepCCEPinternal1{condOfInt}(1:length(tCCEP),:,j));
    tempDbsBase1 = mean(tempDbs1(tCCEP<-0.25,:),1);
    if subtractPre
        tempDbsNormalized1 = tempDbs1 - repmat(tempDbsBase1,[size(tempDbs1,1),1]);
    else
        tempDbsNormalized1 = tempDbs1;
    end
    
    % opposite polarity
    tempDbs2 = squeeze(DBS_sepCCEPinternal2{condOfInt}(1:length(tCCEP),:,j));
    tempDbsBase2 = mean(tempDbs2(tCCEP<-0.25,:),1);
    if subtractPre
        tempDbsNormalized2 = tempDbs2 - repmat(tempDbsBase2,[size(tempDbs2,1),1]);
    else
        tempDbsNormalized2 = tempDbs2;
    end
    
    % normalize if one is longer than the other
    if size(tempDbsNormalized2,2)>size(tempDbsNormalized1,2)
        tempDbsNormalized2 = tempDbsNormalized2(:,1:size(tempDbsNormalized1,2));
    else
        tempDbsNormalized1 = tempDbsNormalized1(:,1:size(tempDbsNormalized2,2));
    end
    
    
    tempDbsNormalized = [tempDbsNormalized1 tempDbsNormalized2];
    % rereference if wanted
    
    
    mu = mean(tempDbsNormalized,2);
    stdError = std(tempDbsNormalized,[],2)/sqrt(size(tempDbsNormalized,2));
    
    
    %mu = mean(squeeze(DBS_sepCCEPinternal{condOfInt}(1:length(tCCEP),:,j)),2);
    %stdError = std(squeeze(DBS_sepCCEPinternal{condOfInt}(1:length(tCCEP),:,j)),[],2)/sqrt(size(squeeze(DBS_sepCCEPinternal{condOfInt}(1:length(tCCEP),:,j)),2));
    
    switch plotType
        case 'errorPlot'
            plot(tCCEP,mu)
            hold on
            plot(tCCEP, mu+stdError, ':');
            hold on;
            plot(tCCEP, mu-stdError, ':');
        case 'doShade'
            plotBTLError(tCCEP,tempDbsNormalized,'CI');
        case 'individual'
            plot(tCCEP,tempDbsNormalized)
            hold on
            plot(mu,'k','linewidth',2)
    end
    %ylim([-6e-6 6e-6])
    xlim([0 5])
    % put a box around the stimulation channels of interest if need be
    if ismember(j,stimChans)
        ax = gca;
        ax.Box = 'on';
        ax.XColor = 'red';
        ax.YColor = 'red';
        ax.LineWidth = 2;
        title(['DBS Channel ',num2str(j)],'color','red');
        
    else
        title(['DBS Channel ',num2str(j)]);
    end
    
end
ylabel('Voltage (V)')
xlabel('time (ms)')
%%

fig1  = figure;
p = numSubplots(numEco);
p(1) = 2
p(2) = 8
fig2 = figure;

for j = 1:numEco
    figure(fig1);
    subplot(p(1),p(2),j)
    % positive polarity
    tempEco1 = squeeze(ECoG_sepCCEPinternal1{condOfInt}(1:length(tCCEP),:,j));
    tempEcoBase1 = mean(tempEco1(tCCEP<-0.25,:),1);
    if subtractPre
        tempEcoNormalized1 = tempEco1 - repmat(tempEcoBase1,[size(tempEco1,1),1]);
    else
        tempEcoNormalized1 = tempEco1;
    end
    % negative polarity
    tempEco2 = squeeze(ECoG_sepCCEPinternal2{condOfInt}(1:length(tCCEP),:,j));
    tempEcoBase2 = mean(tempEco2(tCCEP<-0.25,:),1);
    if subtractPre
        tempEcoNormalized2 = tempEco2 - repmat(tempEcoBase2,[size(tempEco2,1),1]);
    else
        tempEcoNormalized2 = tempEco2;
    end
    
    % normalize if one is longer than the other
    if size(tempEcoNormalized2,2)>size(tempEcoNormalized1,2)
        tempEcoNormalized2 = tempEcoNormalized2(:,1:size(tempEcoNormalized1,2));
    else
        tempEcoNormalized1 = tempEcoNormalized1(:,1:size(tempEcoNormalized2,2));
    end
    
    % concatenate matrices
    tempEcoNormalized = [tempEcoNormalized1 tempEcoNormalized2];
    
    mu = mean(tempEcoNormalized,2);
    stdError = std(tempEcoNormalized,[],2)/sqrt(size(tempEcoNormalized,2));
    
    %mu = mean(squeeze(ECoG_sepCCEPinternal{condOfInt}(1:length(tCCEP),:,j)),2);
    % stdError = std(squeeze(ECoG_sepCCEPinternal{condOfInt}(1:length(tCCEP),:,j)),[],2)/sqrt(size(squeeze(ECoG_sepCCEPinternal{condOfInt}(1:length(tCCEP),:,j)),2));
    switch plotType
        case 'errorPlot'
            plot(tCCEP,mu)
            hold on
            plot(tCCEP, mu+stdError, ':');
            hold on;
            plot(tCCEP, mu-stdError, ':');
        case 'doShade'
            plotBTLError(tCCEP,tempEcoNormalized1,'CI','r');
            hold on
            plotBTLError(tCCEP,tempEcoNormalized2,'CI');
        case 'individual'
            plot(tCCEP,tempEcoNormalized1,'r')
            hold on
            plot(tCCEP,tempEcoNormalized2)
            plot(mu,'k','linewidth',2)
    end
    
    %ylim([-1e-5 1e-5])
    %xlim([min(tCCEP) 5])
        xlim([-1 3])

    title(['ECoG Channel ',num2str(j)])
    
end
ylabel('Voltage (V)')
xlabel('time (ms)')
p = numSubplots(numDBS);
p(1) = 2
p(2) = 4
for j = 1:numDBS
    figure(fig2);
    subplot(p(1),p(2),j)
    % one polarity
    tempDbs1 = squeeze(DBS_sepCCEPinternal1{condOfInt}(1:length(tCCEP),:,j));
    tempDbsBase1 = mean(tempDbs1(tCCEP<-0.25,:),1);
    if subtractPre
        tempDbsNormalized1 = tempDbs1 - repmat(tempDbsBase1,[size(tempDbs1,1),1]);
    else
        tempDbsNormalized1 = tempDbs1;
    end
    % opposite polarity
    tempDbs2 = squeeze(DBS_sepCCEPinternal2{condOfInt}(1:length(tCCEP),:,j));
    tempDbsBase2 = mean(tempDbs2(tCCEP<-0.25,:),1);
    if subtractPre
        tempDbsNormalized2 = tempDbs2 - repmat(tempDbsBase2,[size(tempDbs2,1),1]);
    else
        tempDbsNormalized2 = tempDbs2;
    end
    
    % normalize if one is longer than the other
    if size(tempDbsNormalized2,2)>size(tempDbsNormalized1,2)
        tempDbsNormalized2 = tempDbsNormalized2(:,1:size(tempDbsNormalized1,2));
    else
        tempDbsNormalized1 = tempDbsNormalized1(:,1:size(tempDbsNormalized2,2));
    end
    
    
    % concatenate matrices
    
    tempDbsNormalized = [tempDbsNormalized1 tempDbsNormalized2];
    mu = mean(tempDbsNormalized,2);
    stdError = std(tempDbsNormalized,[],2)/sqrt(size(tempDbsNormalized,2));
    
    
    %mu = mean(squeeze(DBS_sepCCEPinternal{condOfInt}(1:length(tCCEP),:,j)),2);
    %stdError = std(squeeze(DBS_sepCCEPinternal{condOfInt}(1:length(tCCEP),:,j)),[],2)/sqrt(size(squeeze(DBS_sepCCEPinternal{condOfInt}(1:length(tCCEP),:,j)),2));
    
    switch plotType
        case 'errorPlot'
            plot(tCCEP,mu)
            hold on
            plot(tCCEP, mu+stdError, ':');
            hold on;
            plot(tCCEP, mu-stdError, ':');
        case 'doShade'
            plotBTLError(tCCEP,tempDbsNormalized1,'CI','r');
            hold on
            plotBTLError(tCCEP,tempDbsNormalized2,'CI');
        case 'individual'
            plot(tCCEP,tempDbsNormalized1,'r')
            hold on
            plot(tCCEP,tempDbsNormalized2)
            plot(mu,'k','linewidth',2)
    end
            xlim([-1 3])
    %ylim([-10e-6 10e-6])

    % put a box around the stimulation channels of interest if need be
    if ismember(j,stimChans)
        ax = gca;
        ax.Box = 'on';
        ax.XColor = 'red';
        ax.YColor = 'red';
        ax.LineWidth = 2;
        title(['DBS Channel ',num2str(j)],'color','red');
        
    else
        title(['DBS Channel ',num2str(j)]);
    end
    
end
ylabel('Voltage (V)')
xlabel('time (ms)')
