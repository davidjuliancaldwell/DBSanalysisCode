% Visualize peak-to-peak extraction for specified subjects
% Replicates the paper pipeline: avg 5 trials -> SG smooth -> peak_to_peak
% Shows extracted waveforms with peak/trough markers at each stim level
%
% Configure subjects to run by editing the subjects cell array below.

close all; clear all; clc
Z_ConstantsDBS_PairedPulse;
matlab_dir = 'MATLAB_Converted';
experiment = 'EP_Measurement';

% Match paper pipeline settings
avgTrials = 1;
numAvg = 5;
savePlot = 0;
screenBadChans = 0;
plotCondAvg = 0;
plotPkTr = 0;
tryArtifact = 0;
chanReref = 1;
rerefMode = 'none';
smooth = 1;
order = 3;
framelen = 91;

% --- Configure which subjects to run ---
subjects = {'41a73', '68574'};

for s = 1:length(subjects)
    sid = subjects{s};

    % Subject config: blocks/channels from R config, labels from master script
    switch sid
        case '46c2a'
            blocks = [1 3];
            chanIntList = [6];
            blockNames = {'baseline','post 25 ms A/B'};
        case 'c963f'
            blocks = [2 4 6 8];
            chanIntList = [6];
            blockNames = {'baseline 1','post A/B 200 ms','baseline 2','post A/A 200 ms'};
        case '2e114'
            blocks = [1 3 5 7 10 12];
            chanIntList = [4];
            blockNames = {'baseline 1','post 25 ms A/B','baseline 2',...
                'post 50 ms A/B','baseline 3','post A only'};
        case '3d413'
            blocks = [1 2 3 4 5 6 7 8 9 10];
            chanIntList = [4 6];
            blockNames = {'asleep 1 (ch6)','asleep 2 (ch4)','awake 1 (ch6)',...
                'awake 2 (ch4)','asleep 3 (ch6)','asleep 4 (ch4)',...
                'awake 3 (ch6)','awake 4 (ch4)','asleep 5 (ch6)','asleep 6 (ch4)'};
        case 'fe7df'
            blocks = [1 2 3 4 5];
            chanIntList = [4 5];
            blockNames = {'baseline 1','baseline 2','post A/B 200 ms 1',...
                'post A/B 200 ms 2','post A/B 200 ms 3'};
        case 'e6f3c'
            blocks = [1 2 3 4 5];
            chanIntList = [4 5];
            blockNames = {'baseline 1','baseline 2','post A/B 200 ms 1',...
                'post A/B 200 ms 2','post A/B 200 ms 3'};
        case '9f852'
            blocks = [1 2 3 4 5 6];
            chanIntList = [5 6];
            blockNames = {'baseline 1','post A/B 200 ms 1','baseline 2',...
                'post A/B 25 ms 1','baseline 3','post A/A 200 ms'};
        case '8e907'
            blocks = [1 2 3 4 5 6 7 8];
            chanIntList = [6];
            blockNames = {'baseline 1','baseline 2','post A/B 200 ms',...
                'baseline 3','post A/A 200 ms','baseline 4',...
                'post A/B 25 ms','baseline 5'};
        case '08b13'
            blocks = [1 3 5 6 7 8];
            chanIntList = [5 6];
            blockNames = {'baseline 1','baseline 2','post A/B 200 ms 1',...
                'post A/B 200 ms 2/pre A/A','post A/A 200 ms 1',...
                'post A/A 200 ms 2'};
        case 'e9c9b'
            blocks = [1 2 3 4 5 6 7 8 9 10];
            chanIntList = [5 6];
            blockNames = {'baseline 1','baseline 2','post A/A 200 ms 1',...
                'baseline 3 (post A/A)','post A/B 200 ms 1',...
                'baseline 4 (post A/B)','post A/B 25 ms 1',...
                'baseline 5 (post A/B)','baseline 6','post A/A 25 ms 1'};
        case '41a73'
            blocks = [1 2 3 4 5 6 7 8 9 10];
            chanIntList = [5];
            blockNames = {'baseline 1','baseline 2 (pre cond)','post A/B 200 ms',...
                'baseline 3 (post A/B 200)','baseline 4',...
                'post A/A 200 ms','baseline 5 (post A/A)',...
                'baseline 6','baseline 7','baseline 8'};
        case '68574'
            blocks = [1 2 3 4 5 6 7 8 9 10 11];
            chanIntList = [4 7];
            blockNames = {'baseline 1','baseline 2 (pre cond)',...
                'post A/A 100 ms','baseline 3 (post A/A)',...
                'post A/B 100 ms','baseline 4 (post A/B)',...
                'post A/A 200 ms','baseline 5 (post A/A)',...
                'post A/B 200 ms','baseline 6 (post A/B)',...
                'baseline 7 - pre DBS'};
        case '01fee'
            blocks = [1 2 3 4 5 6 7 8 9 10 11];
            chanIntList = [5 8];
            blockNames = {'baseline 1','baseline 2','post A/B 100 ms 1',...
                'baseline 3 (post A/B)','post A/A 100 ms 1',...
                'baseline 4 (post A/A)','post A/B 200 ms 1',...
                'baseline 5 (post A/B, pre DBS)','during DBS',...
                'post DBS 1','post DBS 2'};
        case 'a23ed'
            blocks = [1 2 3 4 5 6 7 8 9];
            chanIntList = [5 8];
            blockNames = {'baseline 1','baseline 2','post A/B 200 ms (5 min)',...
                'baseline 3 (post A/B)','baseline 4',...
                'post A/B 200 ms (15 min)','baseline 5',...
                'post A/A 200 ms (15 min)','baseline 6'};
        otherwise
            error('Unknown subject: %s', sid);
    end

    fprintf('\n=== Processing %s ===\n', sid);

    % Temporarily set legendText for prepare_EP_blocks compatibility
    legendText = blockNames;

    % Load data
    prepare_EP_blocks

    fprintf('Stim levels: ');
    fprintf('%d ', stimLevelUniq);
    fprintf('\n');

    nStimLevels = length(stimLevelUniq);
    stimLabels = arrayfun(@(x) sprintf('%.1f mA', x/1e3), stimLevelUniq, 'UniformOutput', false);
    cmap = lines(nStimLevels);

    for chanIdx = 1:length(chanIntList)
        chanInt = chanIntList(chanIdx);

        % --- Per-block figures: all stim levels in one 2x2 subplot ---
        for bi = 1:length(blocks)
            figure('Units','inches','Position',[1 1 12 8]);
            sgtitle(sprintf('%s - Block %d: %s, Chan %d (avg %d trials)', ...
                sid, blocks(bi), blockNames{bi}, chanInt, numAvg), 'fontsize', 14);

            for condInt = 1:nStimLevels
                epochData = epochsEPblock{bi}{condInt}; % Time x Channels x Trials
                chanData = squeeze(epochData(:, chanInt, :)); % Time x Trials

                % Average every numAvg trials (paper pipeline)
                chanDataAvg = avg_every_p_elems(chanData, numAvg);
                nAvgTrials = size(chanDataAvg, 2);

                % Window to tBegin-tEnd
                tMask = tEpoch > tBegin & tEpoch < tEnd;
                tWindowed = tEpoch(tMask);

                subplot(2, 2, condInt);
                hold on

                ppVals = [];
                for tr = 1:nAvgTrials
                    trialData = chanDataAvg(tMask, tr);

                    % Savitzky-Golay smooth
                    trialSmooth = sgolayfilt_complete(trialData, order, framelen);

                    % Peak to peak extraction (match main pipeline: empty -> NaN)
                    [amp, pk_loc, tr_loc] = peak_to_peak(trialSmooth);
                    if isempty(amp)
                        amp = nan; pk_loc = nan; tr_loc = nan;
                    end

                    % Plot smoothed waveform
                    plot(tWindowed, 1e6*trialSmooth, 'color', [0.6 0.6 0.6], 'linewidth', 0.5);

                    % Mark peak and trough
                    if ~isnan(amp)
                        plot(tWindowed(pk_loc), 1e6*trialSmooth(pk_loc), 'v', ...
                            'color', cmap(condInt,:), 'markersize', 6, 'markerfacecolor', cmap(condInt,:));
                        plot(tWindowed(tr_loc), 1e6*trialSmooth(tr_loc), '^', ...
                            'color', cmap(condInt,:), 'markersize', 6, 'markerfacecolor', cmap(condInt,:));
                        ppVals(end+1) = amp * 1e6;
                    end
                end

                % Plot the grand mean (of all averaged trials)
                meanSmooth = sgolayfilt_complete(mean(chanDataAvg(tMask,:), 2), order, framelen);
                plot(tWindowed, 1e6*meanSmooth, 'color', cmap(condInt,:), 'linewidth', 2.5);

                xlim([tBegin-1 tEnd+2])
                ylim([-350 350])
                ylabel('Voltage (\muV)')
                xlabel('time (ms)')
                if ~isempty(ppVals)
                    title(sprintf('%s  |  PP: %.0f \\pm %.0f \\muV (n=%d)', ...
                        stimLabels{condInt}, median(ppVals), std(ppVals), length(ppVals)));
                else
                    title(sprintf('%s  |  no valid PP', stimLabels{condInt}));
                end
                set(gca, 'fontsize', 10)
            end
        end

        % --- Summary boxplot: PP values across blocks and stim levels ---
        figure('Units','inches','Position',[1 1 10 6]);

        allPP = [];
        allStim = [];

        for bi = 1:length(blocks)
            for condInt = 1:nStimLevels
                epochData = epochsEPblock{bi}{condInt};
                chanData = squeeze(epochData(:, chanInt, :));
                chanDataAvg = avg_every_p_elems(chanData, numAvg);
                tMask = tEpoch > tBegin & tEpoch < tEnd;

                for tr = 1:size(chanDataAvg, 2)
                    trialSmooth = sgolayfilt_complete(chanDataAvg(tMask, tr), order, framelen);
                    [amp, ~, ~] = peak_to_peak(trialSmooth);
                    if isempty(amp), amp = nan; end
                    if ~isnan(amp)
                        allPP(end+1) = amp * 1e6;
                        allStim(end+1) = condInt;
                    end
                end
            end
        end

        boxplot(allPP, allStim, 'Labels', stimLabels);
        ylabel('Peak-to-Peak Amplitude (\muV)')
        xlabel('Stimulation Level')
        title(sprintf('%s Chan %d - PP amplitude by stim level (avg %d, all R blocks)', sid, chanInt, numAvg))
        set(gca, 'fontsize', 12)
    end
end
