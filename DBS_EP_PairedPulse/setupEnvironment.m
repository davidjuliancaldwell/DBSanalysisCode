% setupEnvironment.m — Set environment variables and paths for DBS paired pulse analysis.
% Run this once at the start of a MATLAB session before calling any analysis scripts.
% Modify the paths below to match your machine.

if ispc
    setenv('dbs_subject_dir', 'G:\My Drive\GRIDLabDavidShared\DBS\');
    setenv('OUTPUT_DIR', 'C:\Users\djcald.CSENETID\Data\Output');

elseif ismac
    setenv('dbs_subject_dir', '/Users/davidcaldwell/Library/CloudStorage/OneDrive-UCSF/UWGoogleDriveBackup/djcald_uw_backup_v2/GRIDLabDavidShared_v2/DBS');
    setenv('OUTPUT_DIR', fullfile(pwd, 'output_temp'));
end

% Create OUTPUT_DIR subdirectory if needed
outSubdir = fullfile(getenv('OUTPUT_DIR'), 'DBS', 'pairedPulse');
if ~exist(outSubdir, 'dir')
    mkdir(outSubdir);
end

% Add this repo and MATLAB_ECoG_code to the path
addpath(genpath(fileparts(pwd)));
addpath(genpath('/Users/davidcaldwell/code/MATLAB_ECoG_code'));

fprintf('Environment configured.\n');
fprintf('  dbs_subject_dir = %s\n', getenv('dbs_subject_dir'));
fprintf('  OUTPUT_DIR      = %s\n', getenv('OUTPUT_DIR'));
