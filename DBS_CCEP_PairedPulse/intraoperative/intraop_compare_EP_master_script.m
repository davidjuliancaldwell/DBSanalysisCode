%% script to look at EPs intraoperatively after each session is complete
%
% David.J.Caldwell 10.10.2018
% script to analyze EP results from paired pulse DBS experiments

close all;clear all;clc % prepare workspace

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
needToConvert = 1; % if needing to run conversion script

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if needToConvert
    ConvertTDTRecordingToMAT_v2_loop
end

%% set variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% here are subjects that we have acquired MEP data on
fileFolder = 'C:\Subjects';% file folder
sid = 'e6f3c';
matlabFolder = 'MATLAB_Converted';
dataFolder = [sid '_EP_Measurement'];

% when to start looking for peak to peak values
tBegin = 2.1; % ms
tEnd = 35; % ms %

blockCount = 1; % where to start counting (should be 1)
blocks = [1 4]; % which of the converted blocks to do
stimChansVec = [5 6; 5 6]; % stimChans for each block

plotCondAvg = 0; % plot average of each condition as we go, true/false
savePlot = 0; % true/false save plot
legendText = {'pre conditioning' ,'post 25 ms conditioning 1'}; % what is the legend text

fprintf([sid,'\n'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% prepare data
prepare_EP_blocks_intraop

%% compare multiples blocks

analyze_EP_compare_multiple_blocks_intraop