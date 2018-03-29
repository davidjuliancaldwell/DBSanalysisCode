% Script to compare conditions per subject
close all;clear all;clc
Z_ConstantsDBS_externalEPs

sid = 'b26b7';
file = 1;
plotIt = 1;
switch sid
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

plot_EPs(ucondition,DBS_sep,ECoG_sep,t,stimChans,windowInt)