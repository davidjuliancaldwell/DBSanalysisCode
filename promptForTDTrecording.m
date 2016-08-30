function [dataStructure,filepath] = promptForTDTrecording(baseDir)
% function filepath = promptForTDTrecording(baseDir)
%
% opens a file selection dialog specific to any .mat data file (primarily
% converted from the TDT)
    if(~exist('baseDir', 'var'))        
        baseDir = myGetenv('subject_dir');
    end
    
    currentDir = pwd;

    try
        cd(baseDir);
        [name,path] = uigetfile('*.mat','MultiSelect', 'off');
    catch
        cd(currentDir);
        filepath = '';
        return;
    end

    cd(currentDir);

    filepath = fullfile(path,name);
    dataStructure = load(filepath);
end
