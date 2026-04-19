function TouchDir(path)
% Create the directory at `path`, plus every intermediate directory,
% idempotently. Copied from MATLAB_ECoG_code/DataPrep/TouchDir.m so
% DBSanalysisCode's local SaveFig has a local dependency.

slashies = find(path == '/' | path == '\');
if isempty(slashies) || slashies(end) ~= length(path)
    slashies(end+1) = length(path);
end

for i = slashies
    destDir = path(1:i);
    [~, ~, ~] = mkdir(destDir);
end
