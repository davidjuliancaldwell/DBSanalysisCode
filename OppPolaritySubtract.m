%% DJC 11-2-2016 - updated 6-6-2017 
% Script to subtract opposite stimulation waveform pairs

Z_ConstantsDBS

prompt = {'What is the condition of interest?'};
dlg_title = 'Condition of interest ';
num_lines = 1;
defaultans = {'4'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

condOfInt = str2num(answer{1});

load(fullfile(OUTPUT_DIR,'stimInternal_l_singleDBS_3_4_fs_185.mat'));
ECoG_sepCCEPinternal1 = ECoG_sepCCEPinternal;
DBS_sepCCEPinternal1 = DBS_sepCCEPinternal;


load(fullfile(OUTPUT_DIR,'stimInternal_l_singleDBS_4_3_fs_185.mat'));
ECoG_sepCCEPinternal2 = ECoG_sepCCEPinternal;
DBS_sepCCEPinternal2 = DBS_sepCCEPinternal;

figure

numEco = size(ECoG_sepCCEPinternal{1},3);
numDBS = size(DBS_sepCCEPinternal{1},3);
%%
for j = 1:numEco
    subplot(4,4,j)
    % positive polarity
    tempEco1 = squeeze(ECoG_sepCCEPinternal1{condOfInt}(1:length(tCCEP),:,j));
    tempEcoBase1 = mean(tempEco1(tCCEP<-0.25,:),1);
    tempEcoNormalized1 = tempEco1;
    %tempEcoNormalized1 = tempEco1 - repmat(tempEcoBase1,[size(tempEco1,1),1]);
    
    % negative polarity
    tempEco2 = squeeze(ECoG_sepCCEPinternal2{condOfInt}(1:length(tCCEP),:,j));
    tempEcoBase2 = mean(tempEco2(tCCEP<-0.25,:),1);
    tempEcoNormalized2 = tempEco2;
    %tempEcoNormalized2 = tempEco2 - repmat(tempEcoBase2,[size(tempEco2,1),1]);
    
    
    tempEcoNormalized = [tempEcoNormalized1 tempEcoNormalized2];

    
    mu = mean(tempEcoNormalized,2);
    stdError = std(tempEcoNormalized,[],2)/sqrt(size(tempEcoNormalized,2));
    
    %mu = mean(squeeze(ECoG_sepCCEPinternal{condOfInt}(1:length(tCCEP),:,j)),2);
    % stdError = std(squeeze(ECoG_sepCCEPinternal{condOfInt}(1:length(tCCEP),:,j)),[],2)/sqrt(size(squeeze(ECoG_sepCCEPinternal{condOfInt}(1:length(tCCEP),:,j)),2));
    plot(tCCEP,mu)
    hold on
    plot(tCCEP, mu+stdError, ':');
    hold on;
    plot(tCCEP, mu-stdError, ':');
    ylabel('Voltage (V)')
    xlabel('time (ms)')
    ylim([-1e-6 1e-6])
    xlim([min(tCCEP) 5])
    
    title(['ECoG Channel ',num2str(j)])
    
end

figure
for j = 1:numDBS
    
    subplot(2,4,j)
    % one polarity
    tempDbs1 = squeeze(DBS_sepCCEPinternal1{condOfInt}(1:length(tCCEP),:,j));
    tempDbsBase1 = mean(tempDbs1(tCCEP<-0.25,:),1);
    %tempDbsNormalized1 = tempDbs1 - repmat(tempDbsBase1,[size(tempDbs1,1),1]);
    tempDbsNormalized1 = tempDbs1;
    
    % opposite polarity 
    tempDbs2 = squeeze(DBS_sepCCEPinternal2{condOfInt}(1:length(tCCEP),:,j));
    tempDbsBase2 = mean(tempDbs2(tCCEP<-0.25,:),1);
    %tempDbsNormalized2 = tempDbs2 - repmat(tempDbsBase2,[size(tempDbs2,1),1]);
    tempDbsNormalized2 = tempDbs2;
    
    tempDbsNormalized = [tempDbsNormalized1 tempDbsNormalized2];
    
    mu = mean(tempDbsNormalized,2);
    stdError = std(tempDbsNormalized,[],2)/sqrt(size(tempDbsNormalized,2));
    
    
    %mu = mean(squeeze(DBS_sepCCEPinternal{condOfInt}(1:length(tCCEP),:,j)),2);
    %stdError = std(squeeze(DBS_sepCCEPinternal{condOfInt}(1:length(tCCEP),:,j)),[],2)/sqrt(size(squeeze(DBS_sepCCEPinternal{condOfInt}(1:length(tCCEP),:,j)),2));
    
    plot(tCCEP,mu)
    hold on
    plot(tCCEP, mu+stdError, ':');
    hold on;
    plot(tCCEP, mu-stdError, ':');
    ylabel('Voltage (V)')
    xlabel('time (ms)')
    ylim([-1e-4 1e-4])
    xlim([min(tCCEP) 5])
    
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

