% figure
% data = [data1;data2;data3;data4;data5;data6;data7;data8];
% for i=1:8
%     subplot(4,2,1)
%     plot(data

figure
subplot(4,2,1)
plot(data1), title('1')
subplot(4,2,2)
plot(data2), title('2')
subplot(4,2,3)
plot(data3), title('3')
subplot(4,2,4)
plot(data4), title('4')
subplot(4,2,5)
plot(data5), title('5')
subplot(4,2,6)
plot(data6), title('6')
subplot(4,2,7)
plot(data7), title('7')
subplot(4,2,8)
plot(data8), title('8')

%%
figure
subplot(8,1,1)
plot(data1), title('1')
subplot(8,1,2)
plot(data2), title('2')
subplot(8,1,3)
plot(data3), title('3')
subplot(8,1,4)
plot(data4), title('4')
subplot(8,1,5)
plot(data5), title('5')
subplot(8,1,6)
plot(data6), title('6')
subplot(8,1,7)
plot(data7), title('7')
subplot(8,1,8)
plot(data8), title('8')

%%
 [num_data text_data] = xlsread('C:\Users\jcronin\Google Drive\GRIDLabDavidShared\DBS\be99a\MEP files, EMG for MEP\Brickert_Left Cranium MEP_THENAR.csv',...
    1, 'B8:ATF8');
%%
for i=1:length(text_data)
    trials(i) = regexprep(text_data(i), '/.*', '');
end

trials = str2double(trials);
dx = diff(trials);
ind = find(dx~=0)+1;
amp = trials(ind);

T = [ind', amp'];

%% Plot TDT data
figure
subplot(2,1,1)
plot(Sing.data(:,1))
ylabel('DBS output value (V)')
title('TDT Trial 1')
subplot(2,1,2)
plot(ECOG.data(:,6))
ylabel('ECoG ch. 6 value (V)')

%% Plot DBS electrodes
figure
subplot(4,1,1), plot(DBSs.data(:,1))
subplot(4,1,2), plot(DBSs.data(:,2))
subplot(4,1,3), plot(DBSs.data(:,3))
subplot(4,1,4), plot(DBSs.data(:,4))

%%
figure
subplot(4,1,1), plot(DBSs.data(:,5))
subplot(4,1,2), plot(DBSs.data(:,6))
subplot(4,1,3), plot(DBSs.data(:,7))
subplot(4,1,4), plot(DBSs.data(:,8))