%% script to look at EPs intraoperatively after each session is complete
%
% David.J.Caldwell 10.10.2018
% script to analyze EP results from paired pulse DBS experiments

%close all;clear all;clc % prepare workspace

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
needToConvert = 1; % if needing to run conversion script

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if needToConvert
    ConvertTDTRecordingToMat_v2_loop
end

%% set variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% here are subjects that we have acquired MEP data on
fileFolder = 'C:\Subjects';% file folder
sid = '010dcb';
matlabFolder = 'MATLAB_Converted';
dataFolder = ['EP_Measurement'];

% when to start looking for peak to peak values
tBegin = 5; % ms
tEnd = 60; % ms %

blockCount = 1; % where to start counting (should be 1)
blocks = [1]; % which of the converted blocks to do
stimChansVec = [7 6]; % stimChans for each block
chanIntList = [5];% which channels to look at
plotCondAvg = 0; % plot average of each condition as we go, true/false
savePlot = 0; % true/false save plot
%legendText = {' pre conditioning 1',' pre conditioning 2',...
 %   'post A/A 100 ms conditioning 1', '2nd baseline','post A/B 100',...
%    '3rd baseline','post A/A 200 ms','4th baseline','post A/B 200 ms'};
%legendText = {'baseline 1','baseline 2', 'post-conditioning A/B 100ms',...
%    'baseline 3','post-conditioning A/A 100ms','baseline 4', 'post A/B 200 ms','baseline 5'};
legendText = {'post 200 ms A/B 5 mins','baseline 3','baseline 4','post 200 ms A/B 15 mins'} ;
plotCI = 0; % plot shaded confidence intervals, this is what was crashing
fprintf([sid,'\n'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare data
prepare_EP_blocks_intraop_epilepsy
% compare multiples blocks
%%
analyze_EP_compare_multiple_blocks_intraop_epilepsy