function plot_EPs(ucondition,DBSSep,ECoGSep,t,stimChans,windowInt)
%PLOT_EPS Summary of this function goes here
%   Detailed explanation goes here
% plot ECoG Electrodes


for i = 1:length(ucondition)
    figure
    ECoGTemp = ECoGSep{i};
    numEco = size(ECoGTemp,2);
    for j = 1:numEco
        subplot(4,4,j);
        plot(t(t>=windowInt(1) & t<=windowInt(2)),squeeze(ECoGTemp(t>=windowInt(1) & t<=windowInt(2),j,:)));
        title(['ECoG Channel ' num2str(j)]);
    end
    xlabel('time (ms)')
    ylabel('voltage (V)')
    %subtitle(['ECoG Electrodes, Condition ' num2str(i)]);
    
end

% plot DBS Electrodes

for i = 1:length(ucondition)
    figure
    DBSTemp = DBSSep{i};
    numDBS = size(DBSTemp,2);
    
    for j = 1:numDBS
        subplot(4,2,j);
        plot(t(t>=windowInt(1) & t<=windowInt(2)),squeeze(DBSTemp(t>=windowInt(1) & t<=windowInt(2),j,:)));
        title(['DBS Channel ' num2str(j)]);
    end
    
    xlabel('time (ms)')
    ylabel('voltage (V)')
  %  subtitle(['DBS Electrodes, Condition ' num2str(i)]);
    
end

% plot channel of interest

% ui box for input
prompt = {'ECoG Channel of interest?','DBS Channel of interest','Condition of Interest?'};
dlgTitle = 'Channel of Interest';
numLines = 1;
defaultans = {'8','1','4'};
answer = inputdlg(prompt,dlgTitle,numLines,defaultans);

ecogChanInt = str2num(answer{1});
DBSChanInt = str2num(answer{2});
condInt = str2num(answer{3});

ECoGTemp = ECoGSep{condInt};
DBSTemp = DBSSep{condInt};

figure
plot(t(t>=windowInt(1) & t<=windowInt(2)),squeeze(ECoGTemp(t>=windowInt(1) & t<=windowInt(2),ecogChanInt,:)))

xlabel('time (ms)')
ylabel('voltage (V)')
title(['ECoG Channel ', num2str(ecogChanInt), ' for Condition ', num2str(condInt)]);

figure
plot(t(t>=windowInt(1) & t<=windowInt(2)),squeeze(DBSTemp(t>=windowInt(1) & t<=windowInt(2),DBSChanInt,:)))

xlabel('time (ms)')
ylabel('voltage (V)')
title(['DBS Channel ', num2str(DBSChanInt), ' for Condition ', num2str(condInt)]);

% look at averages for condition of interest
DBSAveCond  = squeeze(mean(DBSTemp,3));
ECoGAveCond = squeeze(mean(ECoGTemp,3));

figure
plot(t(t>=windowInt(1) & t<=windowInt(2)),DBSAveCond(t>=windowInt(1) & t<=windowInt(2),:));
xlabel('time (ms)')
ylabel('voltage (V)')
title(['Average DBS recording across channels for condition ', num2str(condInt)])

figure
plot(t(t>=windowInt(1) & t<=windowInt(2)),ECoGAveCond(t>=windowInt(1) & t<=windowInt(2),:));
xlabel('time (ms)')
ylabel('voltage (V)')
title(['Average ECoG recording across channels for condition ', num2str(condInt)])
% look at condition of interest
figure
for j = 1:numEco
    subplot(4,4,j)
    plot(t(t>=windowInt(1) & t<=windowInt(2)),ECoGAveCond(t>=windowInt(1) & t<=windowInt(2),j))
    
    title(['ECoG Channel ',num2str(j)])
    ylim([-5e-5 5e-5])
end
xlabel('time (ms)')
ylabel('Voltage (V)')
%subtitle(['ECoG EP response outside train for condition = ' num2str(cond_int)])

figure
for j = 1:numDBS
    subplot(2,4,j)
    plot(t(t>=windowInt(1) & t<=windowInt(2)),DBSAveCond(t>=windowInt(1) & t<=windowInt(2),j))
    xlabel('time (ms)')
    ylabel('Voltage (V)')
    
    % put a box around the stimulation channels of interest if need be
    if ismember(j-1,stimChans)
        ax = gca;
        ax.Box = 'on';
        ax.XColor = 'red';
        ax.YColor = 'red';
        ax.LineWidth = 2;
        title(['DBS Channel ',num2str(j-1)],'color','red');
        
    else
        title(['DBS Channel ',num2str(j-1)]);
    end
    
end

%subtitle(['DBS EP responses outside train  for condition = ' num2str(cond_int)])

end

