% Script to compare conditions per subject
%close all;clear all;clc
Z_ConstantsDBS_externalEPs

% subjects with good response -
%1dd75
%63ce7
%80301
%B305e

sid = '80301';
file = 1;
plotIt = 1;
switch sid
    case '80301'
        switch file
            case 1
                load(fullfile(OUTPUT_DIR,[sid '_stim_R_bothDBS_5_4.mat']));
            case 2
                load(fullfile(OUTPUT_DIR,[sid '_stim_R_bothDBS_6_5.mat']));
            case 3
                load(fullfile(OUTPUT_DIR,[sid '_stim_R_bothDBS_7_6.mat']));
            case 4
                load(fullfile(OUTPUT_DIR,[sid '_stim_R_bothDBS_4_5.mat']));
            case 5
                load(fullfile(OUTPUT_DIR,[sid '_stim_R_bothDBS_5_6.mat']));
            case 6
                load(fullfile(OUTPUT_DIR,[sid '_stim_R_bothDBS_6_7.mat']));
        end
    case '50ad9'
        switch file
            case 1
                load(fullfile(OUTPUT_DIR,[sid '_stim_L_singleDBS_0_1.mat']));
            case 2
                load(fullfile(OUTPUT_DIR,[sid '_stim_L_singleDBS_1_2.mat']));
            case 3
                load(fullfile(OUTPUT_DIR,[sid '_stim_L_singleDBS_2_3.mat']));
            case 4
                load(fullfile(OUTPUT_DIR,[sid '_stim_L_singleDBS_3_2.mat']));
            case 5
                load(fullfile(OUTPUT_DIR,[sid '_stim_L_singleDBS_2_1.mat']));
            case 6
                load(fullfile(OUTPUT_DIR,[sid '_stim_L_singleDBS_1_0.mat']));
        end
        %     case 'bb908'
        %         switch file
        %             case 1
        %                 load(fullfile(OUTPUT_DIR,[sid '_stim_l_bothDBS_2_3.mat']));
        %             case 2
        %         end
        
    case 'b26b7'
        switch file
            case 1
                load(fullfile(OUTPUT_DIR,[sid '_stim_l_bothDBS_2_3.mat']));
            case 2
                load(fullfile(OUTPUT_DIR,[sid '_stim_l_bothDBS_3_2.mat']));
        end
end
%%

windowInt = [500 1250];
plot_EPs(ucondition,DBS_sep,ECoG_sep,t,stimChans,windowInt);
%%
tBegin = 510;
tEnd = 575;
[ECoGPP,DBSPP,~,~] = extract_PP(ucondition,DBS_sep,ECoG_sep,t,stimChans,tBegin,tEnd);

figure
elecs = 1:size(ECoGPP,1);
for i = ucondition
    
    plot(elecs,ECoGPP(:,i),'o')
    
end
legend({'1','2','3','4'});
xlabel('electrodes')
ylabel('peak to peak')

%% look at fft of DBS electrodes as confirmation
[f,p1] = fft_compute(dbs_fs,DBS_sep{1}(t>500,:,:));

p1Avg = mean(p1,3);

figure
for chan = 1:8
    subplot(8,1,chan)
    plot(f,p1Avg(:,chan))
    xlim([0 200])
    ylim([0 5e-6])
end