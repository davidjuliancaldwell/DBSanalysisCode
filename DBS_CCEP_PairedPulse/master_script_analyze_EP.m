% script to analyze EP results from paired pulse DBS experiments
close all;clear all;clc
Z_ConstantsDBS_PairedPulse;
%% load in subject
% here are subjects that we have acquired MEP data on

matlab_dir = 'MATLAB_Converted';
experiment = 'EP_Measurement';

% select blocks to look at 
%blocks = [6 7 9 11 13 15]; % this is for the 7/8 pair fo re6f3c
%blocks = [5 10 14 17]; % this is for the 7/8 pair fo re6f3c
blocks = [2 3 5 7 10 11];

% counter - dont change 
blockCount = 1;
% select channel of interest
chanIntList = [3 4];

%sid = '3809e';
%sid = '46c2a';
%sid = 'c963f';
%sid = '2e114';
%sid = '3d413';
%sid = 'fe7df';
%sid = 'e6f3c';
sid = '9f852';
savePlot = 1;
plotCondAvg = 0;

fprintf([sid,'\n'])

%% prepare data 
prepare_EP_blocks

%% compare multiples blocks 

%legendText = {'baseline' ,'post A/B 200 ms delay ','post A/A 25 ms control','post A/B 25 ms delay'};
legendText = {'pre conditioning (2nd pre baseline)' ,'post 25 ms A/B - first time', 'post 25 ms A/A','post 200 ms A/B','post 200 ms A/B - 12  minutes later','post 25 ms A/B - second time'}; % what is the legend text

analyze_EP_compare_multiple_blocks