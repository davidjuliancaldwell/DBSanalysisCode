% Wrapper to regenerate _PairedPulseData_avg_5.csv using the existing
% master_script pipeline with numAvg=5.
% Saves to R_data_regenerated/ — does NOT overwrite R_data/.
%
% Run from the DBS_EP_PairedPulse directory:
%   cd /Users/davidcaldwell/code/DBSanalysisCode/DBS_EP_PairedPulse
%   regenerate_avg5_data

close all; clear all; clc

% Set environment variables
setenv('dbs_subject_dir', '/Users/davidcaldwell/Library/CloudStorage/OneDrive-UCSF/UWGoogleDriveBackup/djcald_uw_backup_v2/GRIDLabDavidShared_v2/DBS');
setenv('OUTPUT_DIR', fullfile(pwd, 'output_temp'));
if ~exist(fullfile(pwd, 'output_temp', 'DBS', 'pairedPulse'), 'dir')
    mkdir(fullfile(pwd, 'output_temp', 'DBS', 'pairedPulse'));
end

addpath(genpath(fileparts(pwd)));
addpath(genpath('/Users/davidcaldwell/code/MATLAB_ECoG_code'));

Z_ConstantsDBS_PairedPulse;

matlab_dir = 'MATLAB_Converted';
experiment = 'EP_Measurement';

% Override key parameters
avgTrialsVec = [1]';
numAvg = 5;
savePlot = 0;
saveData = 0;  % we save manually below
screenBadChans = 0;
plotCondAvg = 0;
plotPkTr = 0;
tryArtifact = 0;
chanReref = 1;
rerefMode = 'none';

sidVecIterate = {'46c2a','c963f','2e114','3d413','fe7df','e6f3c',...
    '9f852','8e907','08b13','e9c9b','41a73','68574',...
    '01fee','a23ed'};

outputDataDir = fullfile(pwd, 'R_data_regenerated');
if ~exist(outputDataDir, 'dir')
    mkdir(outputDataDir);
end

for avgTrials = avgTrialsVec'
    for sidCell = sidVecIterate
        sid = sidCell{:};
        fprintf('\n=== Processing %s (avgTrials=%d, numAvg=%d) ===\n', sid, avgTrials, numAvg);

        try
            % This sets blocks, chanIntList, etc per subject
            eval(['sid = ''' sid ''';']);

            % Use the same subject config as master_script
            switch sid
                case '46c2a'
                    blocks = [1 3]; chanIntList = [6];
                    legendText = {'baseline','post 25 ms A/B'};
                case 'c963f'
                    blocks = [2 4 6 8]; chanIntList = [6];
                    legendText = {'baseline','post 25 ms A/B','baseline 2','post 50 ms A/B'};
                case '2e114'
                    blocks = [1 3 5 7 10 12]; chanIntList = [4];
                    legendText = {'baseline 1','post 25 ms A/B','baseline 2','post 50 ms A/B','baseline 3','post A only'};
                case '3d413'
                    blocks = [1 2 3 4 5 6 7 8 9 10]; chanIntList = [4 6];
                    legendText = {"asleep","asleep","awake","awake","awake","awake","asleep","asleep","asleep","asleep"};
                case 'fe7df'
                    blocks = [1 2 3 4 5 6 7 8 9]; chanIntList = [6];
                    legendText = {"baseline 1","baseline 2","baseline 3","post 25 ms A/B 1","baseline 4 post 25 ms A/B 2","post 200 ms A/B 1","baseline 5 post 200 ms A/B 2","baseline 6 post 200 ms A/B 3","baseline 7 post 200 ms A/B 4"};
                case 'e6f3c'
                    blocks = [5 6 7 8 9 10 11 12 13 14 15 16 17]; chanIntList = [7 6];
                    legendText = {"baseline 1 6/5 pair","baseline 1 8/7 pair","baseline 2 6/5 pair","baseline 2 8/7 pair","post A/B 200 ms 1 8/7 pair","post A/B 200 ms 1 6/5 pair","baseline 3 8/7 pair","baseline 3 6/5 pair","post A/A 25 ms 8/7 pair","post A/A 25 ms 6/5 pair","post A/B 25 ms 8/7 pair","post A/B 25 ms 6/5 pair - noisy","post A/B 25 ms 6/ pair - noisy"};
                case '9f852'
                    blocks = [2 3 4 5 6 7 10 11 12]; chanIntList = [4];
                    legendText = {'baseline 2 (pre conditioning)','post A/B 25 ms','baseline 3 (post 25 ms)','post A/A 25 ms','baseline 4 (post 25 ms A/A)','post A/B 200 ms','baseline 5 - post A/B 200 ms 12 minutes later','post A/B 25 ms second time','baseline 6'};
                case '8e907'
                    blocks = [1 2 5 6 7]; chanIntList = [4 5];
                    legendText = {'baseline 1','baseline 2','post A/B 200 ms 1','post A/B 200 ms 2','post A/B 200 ms 3'};
                case '08b13'
                    blocks = [1 2 3 5 6 7]; chanIntList = [6];
                    legendText = {'baseline 1','baseline 2','post A/B 200 ms 1','baseline 3 (post A/B)','post A/A 200 ms 1','post A/A 200 ms 2'};
                case 'e9c9b'
                    blocks = [1 2 3 4 5 6 7 8 9 10]; chanIntList = [6];
                    legendText = {'baseline 1','baseline 2','post A/A 200 ms 1','baseline 3 (post A/A)','post A/B 200 ms 1','baseline 4 (post A/B)','post A/B 25 ms 1','baseline 5 (post A/B)','baseline 6','post A/A 25 ms'};
                case '41a73'
                    blocks = [1 2 3 4 5 6 7 8 9 10]; chanIntList = [5];
                    legendText = {'baseline 1','baseline 2','post A/B 200 ms','baseline 3','baseline 4','post A/A 200 ms','baseline 5','baseline 6','baseline 7','baseline 8'};
                case '68574'
                    blocks = [1 2 3 4 5 6 7 8 9 10 11]; chanIntList = [4 7];
                    legendText = {'baseline 1','baseline 2','post A/A 100 ms 1','baseline 3 (post A/A)','post A/B 100 ms 1','baseline 4 (post A/B)','post A/A 200 ms 1','baseline 5 (post A/A)','post A/B 200 ms 1','baseline 6 (post A/B 200 ms 1)','baseline 7 - pre DBS'};
                case '01fee'
                    blocks = [1 2 3 4 5 6 7 8]; chanIntList = [5 8];
                    legendText = {'baseline 1','baseline 2','post A/B 100 ms 1','baseline 3 (post A/B)','post A/A 100 ms 1','baseline 4 (post A/A)','post A/B 200 ms 1','baseline 5 (post A/B, pre DBS)'};
                case 'a23ed'
                    blocks = [1 2 3 4 5 6 7 8 9]; chanIntList = [5];
                    legendText = {'baseline 1','baseline 2','post A/B 200 ms 1 (5 mins)','baseline 3 (post A/B)','baseline 4','post A/B 200 ms (15 mins)','baseline 5','post A/A 200 ms (15 mins)','baseline 6'};
            end

            % Run the EP preparation and extraction
            prepare_EP_blocks

            % Build output table
            PPvec = [];
            blockVecOut = [];
            stimLevelVecOut = [];
            chanVecOut = [];
            PPfromAvgVec = [];

            for ii = 1:size(signalPPblockST,2)
                for jj = 1:size(signalPPblockST{ii},2)
                    tempPP = signalPPblockST{ii}{jj};
                    tempChannel = repmat([1:8]',[1,size(tempPP,2)]);
                    tempBlock = repmat(blocks(ii),size(tempPP));
                    tempStimLevel = repmat(blockLabel{ii}{jj}',[size(tempPP,1),1]);
                    tempPPfromAvgVec = signalPPblockSTfromAvg{ii}{jj};

                    PPvec = [PPvec; tempPP(:)];
                    blockVecOut = [blockVecOut; tempBlock(:)];
                    stimLevelVecOut = [stimLevelVecOut; tempStimLevel(:)];
                    chanVecOut = [chanVecOut; tempChannel(:)];
                    PPfromAvgVec = [PPfromAvgVec; tempPPfromAvgVec(:)];
                end
            end
            sidVecOut = cellstr(repmat(sid,[size(PPvec),1]));

            T = table(PPvec,blockVecOut,stimLevelVecOut,sidVecOut,chanVecOut,PPfromAvgVec);
            T.Properties.VariableNames = {'PPvec','blockVec','stimLevelVec','sidVec','chanVec','PPfromAvgVec'};

            outfile = fullfile(outputDataDir, [sid '_PairedPulseData_avg_5.csv']);
            writetable(T, outfile, 'Delimiter', ',', 'QuoteStrings', true);
            fprintf('  Saved %s (%d rows)\n', outfile, height(T));

            clearvars signalPPblockST signalPPblockSTfromAvg epochsEPblock
            close all

        catch ME
            fprintf('  ERROR: %s\n', ME.message);
        end
    end
end

fprintf('\nDone. Files in: %s\n', outputDataDir);
