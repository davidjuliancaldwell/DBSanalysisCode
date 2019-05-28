% script to analyze EP results from paired pulse DBS experiments
close all;clear all;clc
Z_ConstantsDBS_PairedPulse;
%% load in subject
% here are subjects that we have acquired MEP data on

matlab_dir = 'MATLAB_Converted';
experiment = 'EP_Measurement';

blocks = [1 4];
blockCount = 1;
%sid = '3809e';
%sid = '46c2a';
%sid = 'c963f';
%sid = '2e114';
%sid = '3d413';
sid = '8e907';
savePlot = 1;
plotCondAvg = 1;

fprintf([sid,'\n'])

%% prepare data 
prepare_EP_blocks

%% compare multiples blocks 

legendText = {'pre conditioning' ,'post 25 ms conditioning 1'};
analyze_EP_compare_multiple_blocks