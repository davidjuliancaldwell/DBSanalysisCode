%% - 3-30-2016 - DJC - response timing file generator - updated 8-19-2016 
% this script will generate two text files. One of these has the sample
% number at which each stimulus train will be delivered. the second has the
% condition which should be read in 

prompt = {'Enter subject name','What is the range of ITI?', 'What is the sample rate of the TDT?','Number Of Trials Per Voltage?','How many stimulation Conditions?','Which file number is this?'};
dlg_title = 'Input';
num_lines = 1;
defaultans = {'DBS','[2.25,2.75]','48828.125','15','4','1'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
sid = answer{1};
ITI = str2num(answer{2});
fs = str2num(answer{3});
numTrials = str2num(answer{4});
stimConds = str2num(answer{5});
fileNum = answer{6};


%% make the timing file 

% add 500 ms to the ITI times to account for the 500 ms pulse trains - then
% 2500 off 

ITIlo = ITI(1)+0.5;
ITIhi = ITI(2)+0.5;

% number of trials for each voltage - 30. 
% so 30 * 6 ?

randTimes = unifrnd(ITIlo,ITIhi,stimConds*numTrials,1);


% here the vector is converted to the sample number where the stimulus
% train should start to be delivered 
sample = 1; % start with sample 1 
pts = [];
for i = 1:length(randTimes)
    sample = floor(sample + randTimes(i)*fs);
    pts = [pts; sample];
end



%% make the conditions file

stimRand = [];

for i = 1:numTrials
   
    stimRand = [stimRand randperm(stimConds)];
    
end


%% write these times to file for stim train delivery

filename = sprintf('%s_stimTrainDelivery_%s.txt',sid,fileNum);
fileID = fopen(filename,'w+');
fprintf(fileID,'%d\r\n',pts);
fclose(fileID);

%% write these times to file for condition 

filename = sprintf('%s_condition_%s.txt',sid,fileNum);
fileID = fopen(filename,'w+');
fprintf(fileID,'%d\r\n',stimRand);
fclose(fileID);