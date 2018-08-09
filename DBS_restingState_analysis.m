%% DBS resting state analysis
%
% script to analyze resting state data from DBS patients with combined ECoG
% and DBS electrodes
%
% David.J.Caldwell 8.8.2018

Z_ConstantsDBS_PairedPulse

sid = 'c963f';
savePlot = 0;

sid


matlab_dir = 'MATLAB_Converted';
experiment = 'resting';

block = 1;

% load in tank
switch sid
    case 'c963f'
        switch block
            % pre first time
            case 1
                load(fullfile(SUB_DIR,sid,matlab_dir,experiment,'restingstate-1.mat'));                
        end
end

%%
% load in data 
ECoG = 4*ECOG.data;
fsECoG = ECOG.info.SamplingRateHz;


DBS = 4*DBSs.data;
fsDBS = DBSs.info.SamplingRateHz;

%%

figure
ax(1) = subplot(2,1,1);
t = 1e3*[0:size(DBS,1)-1]/fsDBS;
plot(t,1e6*highpass(DBS,1,fsDBS))
ylim([-100 100])
ylabel('voltage (\muV)')
xlabel('time (ms)')
title('DBS')


ax(2) = subplot(2,1,2);
plot(t,1e6*highpass(ECoG,1,fsECoG))
ylim([-100 100])
ylabel('voltage (\muV)')
xlabel('time (ms)')
title('ECoG')
linkaxes(ax,'xy')

%%
DBSreref = rereference_CAR_median(DBS,'mean');
ECoGreref = rereference_CAR_median(ECoG,'mean');

figure
ax(1) = subplot(2,1,1);
t = 1e3*[0:size(DBS,1)-1]/fsDBS;
plot(t,1e6*highpass(DBSreref,1,fsDBS))
ylim([-100 100])
ylabel('voltage (\muV)')
xlabel('time (ms)')
title('DBS')


ax(2) = subplot(2,1,2);
plot(t,1e6*highpass(ECoGreref,1,fsECoG))
ylim([-100 100])
ylabel('voltage (\muV)')
xlabel('time (ms)')
title('ECoG')
linkaxes(ax,'xy')




