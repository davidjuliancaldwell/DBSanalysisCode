%% average plots of internal CCEP

prompt = {'What is the condition of interest?'};
dlg_title = 'Condition of interest ';
num_lines = 1;
defaultans = {'4'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

condOfInt = str2num(answer{1});

figure
for j = 1:numEco
    subplot(4,4,j)
    mu = mean(squeeze(ECoG_sepCCEPinternal{condOfInt}(1:length(tCCEP),:,j)),2);
    stdError = std(squeeze(ECoG_sepCCEPinternal{condOfInt}(1:length(tCCEP),:,j)),[],2)/sqrt(size(squeeze(ECoG_sepCCEPinternal{condOfInt}(1:length(tCCEP),:,j)),1));
    plot(tCCEP,mu)
    hold on
    plot(tCCEP, mu+stdError, ':');
    hold on;
    plot(tCCEP, mu-stdError, ':');
    ylabel('Voltage (V)')
    xlabel('time (ms)')
    
    title(['Channel ',num2str(j)])
    
end
subtitle(['ECoG CCEP responses within train for condition = ' num2str(condOfInt)])

figure
for j = 1:numDBS
    subplot(2,4,j)
    mu = mean(squeeze(DBS_sepCCEPinternal{condOfInt}(1:length(tCCEP),:,j)),2);
    stdError = std(squeeze(DBS_sepCCEPinternal{condOfInt}(1:length(tCCEP),:,j)),[],2)/sqrt(size(squeeze(DBS_sepCCEPinternal{condOfInt}(1:length(tCCEP),:,j)),1));
    
    plot(tCCEP,mu)
    hold on
    plot(tCCEP, mu+stdError, ':');
    hold on;
    plot(tCCEP, mu-stdError, ':');
    ylabel('Voltage (V)')
    xlabel('time (ms)')
    
    title(['Channel ',num2str(j)])
end
subtitle(['DBS CCEP responses within train for condition = ' num2str(condOfInt)])

return
