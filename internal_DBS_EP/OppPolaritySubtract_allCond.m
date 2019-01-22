%% Script to subtract opposite stimulation waveform pairs to look at internal CEPS
% Script to subtract opposite stimulation waveform pairs
% David.J.Caldwell 11-2-2016 - updated 6-6-2017
% Updated 1.15.19
close all;clear all;clc
Z_Constants_internal_EP_DBS

sid = '1dd75';
DATA_DIR = fullfile(DATA_DIR_BASE, sid);

load(fullfile(DATA_DIR,'stimInternal_l_singleDBS_3_4_fs_185.mat'));
%load(fullfile(OUTPUT_DIR,'stimInternal_R_bothDBS_5_6_fs_185.mat'));
%load(fullfile(OUTPUT_DIR,'stimInternal_R_singleDBS_1_2_fs_185.mat'));

ECoG_sepCCEPinternal = cellfun(@(x) x*4,ECoG_sepCCEPinternal,'un',0);
DBS_sepCCEPinternal = cellfun(@(x) x*4,DBS_sepCCEPinternal,'un',0);

ECoG_sepCCEPinternal1Raw = ECoG_sepCCEPinternal;
DBS_sepCCEPinternal1Raw = DBS_sepCCEPinternal;

load(fullfile(DATA_DIR,'stimInternal_l_singleDBS_4_3_fs_185.mat'));
%load(fullfile(OUTPUT_DIR,'stimInternal_R_bothDBS_6_5_fs_185.mat'));
%load(fullfile(OUTPUT_DIR,'stimInternal_R_singleDBS_2_1_fs_185.mat'));



ECoG_sepCCEPinternal = cellfun(@(x) x*4,ECoG_sepCCEPinternal,'un',0);
DBS_sepCCEPinternal = cellfun(@(x) x*4,DBS_sepCCEPinternal,'un',0);

ECoG_sepCCEPinternal2Raw = ECoG_sepCCEPinternal;
DBS_sepCCEPinternal2Raw = DBS_sepCCEPinternal;
%%

prompt = {'what are the stim channels? e.g. [1 2]',...
    'what are the channels to exclude from rerereferencing e.g [1 2 3]','"median", "mean","singleChan","bipolar", or "n" rereference for ECoG','"median", "mean", "singleChan","bipolar", or "n" rereference for DBS'};
dlg_title = 'Condition of interest ';
num_lines = 1;
defaultans = {'[1 2]','[9:16]','mean','n'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

stimChans = str2num(answer{1});
badChans = str2num(answer{2});
rrEco = answer{3};
rrDbs = answer{4};

numEco = size(ECoG_sepCCEPinternal{1},3);
numEco = 8;
numDBS = size(DBS_sepCCEPinternal{1},3);
numDBS = 4;
numConds = length(ECoG_sepCCEPinternal);

% cell of conditions
cellConds = {'median','mean','singleChan','bipolar','bipolarPair','none'};
%cellConds = {'bipolar'};
%% rereference

for specificCond = cellConds
    specificCond = specificCond{:};
    rrEco = specificCond;
    rrDbs = 'n';
    
    ecoFig = figure('units','normalized','outerposition',[0 0 1 1]);
   dbsFig = figure('units','normalized','outerposition',[0 0 1 1]);
%        ecoFig = figure('Position', get(0, 'Screensize'));
%     dbsFig = figure('Position', get(0, 'Screensize'));


    for condOfInt = 1:numConds
        %% rereference
        switch rrEco
            case 'median'
                
                %%%%%%%%%%%%%%%%%%%%%%%
                tempEco1 = squeeze(ECoG_sepCCEPinternal1Raw{condOfInt}(1:length(tCCEP),:,:));
                tempEcoNormalized = rereference_CAR_median(tempEco1,rrEco,badChans,[1 3 2]);
                ECoG_sepCCEPinternal1{condOfInt}(1:length(tCCEP),:,:) = tempEcoNormalized;
                %%%%%%%%%%%%%%%%%%%%%%%
                tempEco2 = squeeze(ECoG_sepCCEPinternal2Raw{condOfInt}(1:length(tCCEP),:,:));
                tempEcoNormalized = rereference_CAR_median(tempEco2,rrEco,badChans,[1 3 2]);
                ECoG_sepCCEPinternal2{condOfInt}(1:length(tCCEP),:,:) = tempEcoNormalized;
                
                
            case 'mean'
                
                %%%%%%%%%%%%%%%%%%%%%%%
                tempEco1 = squeeze(ECoG_sepCCEPinternal1Raw{condOfInt}(1:length(tCCEP),:,:));
                tempEcoNormalized = rereference_CAR_median(tempEco1,rrEco,badChans,[1 3 2]);
                ECoG_sepCCEPinternal1{condOfInt}(1:length(tCCEP),:,:) = tempEcoNormalized;
                %%%%%%%%%%%%%%%%%%%%%%%
                tempEco2 = squeeze(ECoG_sepCCEPinternal2Raw{condOfInt}(1:length(tCCEP),:,:));
                tempEcoNormalized = rereference_CAR_median(tempEco2,rrEco,badChans,[1 3 2]);
                ECoG_sepCCEPinternal2{condOfInt}(1:length(tCCEP),:,:) = tempEcoNormalized;
            case 'singleChan'
                %%%%%%%%%%%%%%%%%%%%%%%
                tempEco1 = squeeze(ECoG_sepCCEPinternal1Raw{condOfInt}(1:length(tCCEP),:,:));
                tempEcoNormalized = rereference_CAR_median(tempEco1,rrEco,badChans,[1 3 2],1);% change to 9 for R side strip!! 1 for L side
                ECoG_sepCCEPinternal1{condOfInt}(1:length(tCCEP),:,:) = tempEcoNormalized;
                %%%%%%%%%%%%%%%%%%%%%%%
                tempEco2 = squeeze(ECoG_sepCCEPinternal2Raw{condOfInt}(1:length(tCCEP),:,:));
                tempEcoNormalized = rereference_CAR_median(tempEco2,rrEco,badChans,[1 3 2],1); % change to 9 for R side strip!! 1 for L side 
                ECoG_sepCCEPinternal2{condOfInt}(1:length(tCCEP),:,:) = tempEcoNormalized;
            case 'bipolar'
                %%%%%%%%%%%%%%%%%%%%%%%
                tempEco1 = squeeze(ECoG_sepCCEPinternal1Raw{condOfInt}(1:length(tCCEP),:,:));
                tempEcoNormalized = rereference_CAR_median(tempEco1,rrEco,badChans,[1 3 2]);
                ECoG_sepCCEPinternal1{condOfInt}(1:length(tCCEP),:,:) = tempEcoNormalized;
                %%%%%%%%%%%%%%%%%%%%%%%
                tempEco2 = squeeze(ECoG_sepCCEPinternal2Raw{condOfInt}(1:length(tCCEP),:,:));
                tempEcoNormalized = rereference_CAR_median(tempEco2,rrEco,badChans,[1 3 2]);
                ECoG_sepCCEPinternal2{condOfInt}(1:length(tCCEP),:,:) = tempEcoNormalized;
            case 'bipolarPair'
                %%%%%%%%%%%%%%%%%%%%%%%
                tempEco1 = squeeze(ECoG_sepCCEPinternal1Raw{condOfInt}(1:length(tCCEP),:,:));
                tempEcoNormalized = rereference_CAR_median(tempEco1,rrEco,badChans,[1 3 2]);
                ECoG_sepCCEPinternal1{condOfInt}(1:length(tCCEP),:,:) = tempEcoNormalized;
                %%%%%%%%%%%%%%%%%%%%%%%
                tempEco2 = squeeze(ECoG_sepCCEPinternal2Raw{condOfInt}(1:length(tCCEP),:,:));
                tempEcoNormalized = rereference_CAR_median(tempEco2,rrEco,badChans,[1 3 2]);
                ECoG_sepCCEPinternal2{condOfInt}(1:length(tCCEP),:,:) = tempEcoNormalized;
            otherwise
                ECoG_sepCCEPinternal1 = ECoG_sepCCEPinternal1Raw;
                ECoG_sepCCEPinternal2 = ECoG_sepCCEPinternal2Raw;
                
        end
        
        switch rrDbs
            case 'median'
                
                
                %%%%%%%%%%%%%%%%%%%%%%%
                tempDbs1 = squeeze(DBS_sepCCEPinternal1Raw{condOfInt}(1:length(tCCEP),:,:));
                tempDbsNormalized = rereference_CAR_median(tempDbs1,rrDbs,stimChans,[1 3 2]);
                DBS_sepCCEPinternal1{condOfInt}(1:length(tCCEP),:,:) = tempDbsNormalized;
                %%%%%%%%%%%%%%%%%%%%%%%
                tempDbs2 = squeeze(DBS_sepCCEPinternal2Raw{condOfInt}(1:length(tCCEP),:,:));
                tempDbsNormalized = rereference_CAR_median(tempDbs2,rrDbs,stimChans,[1 3 2]);
                DBS_sepCCEPinternal2{condOfInt}(1:length(tCCEP),:,:) = tempDbsNormalized;
                
            case 'mean'
                
                %%%%%%%%%%%%%%%%%%%%%%%
                tempDbs1 = squeeze(DBS_sepCCEPinternal1Raw{condOfInt}(1:length(tCCEP),:,:));
                tempDbsNormalized = rereference_CAR_median(tempDbs1,rrDbs,stimChans,[1 3 2]);
                DBS_sepCCEPinternal1{condOfInt}(1:length(tCCEP),:,:) = tempDbsNormalized;
                %%%%%%%%%%%%%%%%%%%%%%%
                tempDbs2 = squeeze(DBS_sepCCEPinternal2Raw{condOfInt}(1:length(tCCEP),:,:));
                tempDbsNormalized = rereference_CAR_median(tempDbs2,rrDbs,stimChans,[1 3 2]);
                DBS_sepCCEPinternal2{condOfInt}(1:length(tCCEP),:,:) = tempDbsNormalized;
            case 'singleChan'
                %%%%%%%%%%%%%%%%%%%%%%%
                tempDbs1 = squeeze(DBS_sepCCEPinternal1Raw{condOfInt}(1:length(tCCEP),:,:));
                tempDbsNormalized = rereference_CAR_median(tempDbs1,rrDbs,stimChans,[1 3 2],1);
                DBS_sepCCEPinternal1{condOfInt}(1:length(tCCEP),:,:) = tempDbsNormalized;
                %%%%%%%%%%%%%%%%%%%%%%%
                tempDbs2 = squeeze(DBS_sepCCEPinternal2Raw{condOfInt}(1:length(tCCEP),:,:));
                tempDbsNormalized = rereference_CAR_median(tempDbs2,rrDbs,stimChans,[1 3 2],1);
                DBS_sepCCEPinternal2{condOfInt}(1:length(tCCEP),:,:) = tempDbsNormalized;
            case 'bipolar'
                %%%%%%%%%%%%%%%%%%%%%%%
                tempDbs1 = squeeze(DBS_sepCCEPinternal1Raw{condOfInt}(1:length(tCCEP),:,:));
                tempDbsNormalized = rereference_CAR_median(tempDbs1,rrDbs,stimChans,[1 3 2]);
                DBS_sepCCEPinternal1{condOfInt}(1:length(tCCEP),:,:) = tempDbsNormalized;
                %%%%%%%%%%%%%%%%%%%%%%%
                tempDbs2 = squeeze(DBS_sepCCEPinternal2Raw{condOfInt}(1:length(tCCEP),:,:));
                tempDbsNormalized = rereference_CAR_median(tempDbs2,rrDbs,stimChans,[1 3 2]);
                DBS_sepCCEPinternal2{condOfInt}(1:length(tCCEP),:,:) = tempDbsNormalized;
            case 'bipolarPair'
                %%%%%%%%%%%%%%%%%%%%%%%
                tempDbs1 = squeeze(DBS_sepCCEPinternal1Raw{condOfInt}(1:length(tCCEP),:,:));
                tempDbsNormalized = rereference_CAR_median(tempDbs1,rrDbs,stimChans,[1 3 2]);
                DBS_sepCCEPinternal1{condOfInt}(1:length(tCCEP),:,:) = tempDbsNormalized;
                %%%%%%%%%%%%%%%%%%%%%%%
                tempDbs2 = squeeze(DBS_sepCCEPinternal2Raw{condOfInt}(1:length(tCCEP),:,:));
                tempDbsNormalized = rereference_CAR_median(tempDbs2,rrDbs,stimChans,[1 3 2]);
                DBS_sepCCEPinternal2{condOfInt}(1:length(tCCEP),:,:) = tempDbsNormalized;
            otherwise
                DBS_sepCCEPinternal1 = DBS_sepCCEPinternal1Raw;
                DBS_sepCCEPinternal2 = DBS_sepCCEPinternal2Raw;
                
        end
        %%
        % line vs. shaded
        
        plotTypeVec = {'doShade','error','individual'};
        plotType = plotTypeVec{1};
        
        subtractPre = 0;
        
        figure(ecoFig)
        
        for j = 1:numEco
            subplot(numConds,numEco,j+((condOfInt-1)*numEco))
            %j = j+8; % this is to get right side!
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
            
            ylim([-12e-6 12e-6])
            xlim([0 5])
            
            %xlim([min(tCCEP) 5])
            
            title(['ECoG Channel ',num2str(j)])
            if mod(j,numEco)==1
                ylabel(['Condition ' num2str(condOfInt)])
                
            end
        end
        ylabel('Voltage (V)')
        xlabel('time (ms)')
        figure(dbsFig);
        
        for j = 1:numDBS
            subplot(numConds,numDBS,j+((condOfInt-1)*numDBS))
            %j = j + 4;
            
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
            ylim([-20e-4 20e-4])
            %ylim([-6e-6 6e-6])
            
            xlim([0 5])
            % put a box around the stimulation channels of interest if need be
            if ismember(j,stimChans)
                ax = gca;
                ax.Box = 'on';
                ax.XColor = 'red';
                ax.YColor = 'red';
                ax.LineWidth = 2;
                title(['DBS Channel ',num2str(j-1)],'color','red');
                
            else
                title(['DBS Channel ',num2str(j-1)]);
            end
            
            if mod(j,numDBS)==1
                ylabel(['Condition ' num2str(condOfInt)])
                
            end
            
        end
        ylabel('Voltage (V)')
        xlabel('time (ms)')
        %%
    end
    figure(ecoFig)

    SaveFig(OUTPUT_DIR, sprintf(['walkerECoG-%s-%s-%d-%d'], sid, specificCond, stimChans(1)-1,stimChans(2)-1), 'eps', '-r600');
    SaveFig(OUTPUT_DIR, sprintf(['walkerECoG-%s-%s-%d-%d'], sid, specificCond, stimChans(1)-1,stimChans(2)-1), 'png', '-r600');
    
    figure(dbsFig)

    SaveFig(OUTPUT_DIR, sprintf(['walkerDBS-%s-%s-%d-%d'], sid, specificCond, stimChans(1)-1,stimChans(2)-1), 'eps', '-r600');
    SaveFig(OUTPUT_DIR, sprintf(['walkerDBS-%s-%s-%d-%d'], sid, specificCond, stimChans(1)-1,stimChans(2)-1), 'png', '-r600');
    close all
    clearvars -except ECoG_sepCCEPinternal1Raw ECoG_sepCCEPinternal2Raw DBS_sepCCEPinternal1Raw DBS_sepCCEPinternal2Raw ...
        META_DIR OUTPUT_DIR SIDS numConds numEco numDBS stimChans badChans sid tCCEP
end
return

%%

ecoFigSub = figure('units','normalized','outerposition',[0 0 1 1]);
dbsFigSub = figure('units','normalized','outerposition',[0 0 1 1]);
plotTypeVec = {'doShade','error','individual'};
plotType = plotTypeVec{1};

subtractPre = 0;

for condOfInt = 1:numConds
    
    
    for j = 1:numEco
        figure(ecoFigSub)
        subplot(numConds,numEco,j+((condOfInt-1)*numEco))
        %j = j+8;
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
        
        if mod(j,numEco)==1
            ylabel(['Condition ' num2str(condOfInt)])
            
        end
        
    end
    ylabel('Voltage (V)')
    xlabel('time (ms)')
    
    
    for j = 1:numDBS
        figure(dbsFigSub)
        subplot(numConds,numDBS,j+((condOfInt-1)*numDBS))
            %    j = j+4;

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
        
        
        if mod(j,numDBS)==1
            ylabel(['Condition ' num2str(condOfInt)])
            
        end
        
        
    end
    ylabel('Voltage (V)')
    xlabel('time (ms)')
end
