% this script collects all of the records in a tdt data tank and converts
% it to a mat file
if (exist('myGetenv', 'file'))
    start = myGetenv('subject_dir');    
    if (isempty(start))
        start = pwd;
    end
else
    start = pwd;
end

rawpath = uigetdir(start, 'select a TDT data BLOCK');
outpath = uigetdir(tankpath, 'select an output directory for MAT file');

[tankpath, blockname,ext] = fileparts(rawpath);
experiment = split(blockname,'-');
rawpath_gen = split(rawpath,'-');

num_trials = input('number of trials \n');

for i = 1:num_trials
    
    blocknameOut = strcat(experiment{1},'-',num2str(i));
    %tankpathOut = strcat(rawpath_gen{1},'-',num2str(i));
    mTDT2MAT(tankpath, blocknameOut, outpath);

    
end






