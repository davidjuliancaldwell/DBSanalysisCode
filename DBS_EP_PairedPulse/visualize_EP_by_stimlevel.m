% Visualize mean evoked potentials by stimulation level for 41a73 and 68574
% Uses only blocks included in the R analysis (blockIntLM from subj_ config).
% Plot 1: Mean EP pooled across all R-included blocks, one line per stim level
% Plot 2: Same with CI shading
% Plot 3: One figure per block, all 4 stim levels overlaid

close all; clear all; clc
Z_ConstantsDBS_PairedPulse;
matlab_dir = 'MATLAB_Converted';
experiment = 'EP_Measurement';

avgTrials = 0;
savePlot = 0;
screenBadChans = 0;
plotCondAvg = 0;
plotPkTr = 0;
tryArtifact = 0;
chanReref = 1;
rerefMode = 'none';
numAvg = 3;

subjects = {'41a73', '68574'};

for s = 1:length(subjects)
    sid = subjects{s};

    % Use blocks matching R analysis (blockIntLM from subj_ config files)
    switch sid
        case '41a73'
            blocks = [1 2 3 4 5 6 7 8 9 10];
            chanIntList = [5];
            legendText = {'baseline 1','baseline 2 (pre cond)','post A/B 200 ms',...
                'baseline 3 (post A/B 200)','baseline 4',...
                'post A/A 200 ms','baseline 5 (post A/A)',...
                'baseline 6','baseline 7','baseline 8'};

        case '68574'
            blocks = [1 2 3 4 5 6 7 8 9 10 11];
            chanIntList = [4 7];
            legendText = {'baseline 1','baseline 2 (pre cond)',...
                'post A/A 100 ms','baseline 3 (post A/A)',...
                'post A/B 100 ms','baseline 4 (post A/B)',...
                'post A/A 200 ms','baseline 5 (post A/A)',...
                'post A/B 200 ms','baseline 6 (post A/B)',...
                'baseline 7 - pre DBS'};
    end

    fprintf('\n=== Processing %s ===\n', sid);

    % Run prepare_EP_blocks to load data
    prepare_EP_blocks

    fprintf('Stim levels: ');
    fprintf('%d ', stimLevelUniq);
    fprintf('\n');

    nStimLevels = length(stimLevelUniq);
    stimLabels = arrayfun(@(x) sprintf('%.1f mA', x/1e3), stimLevelUniq, 'UniformOutput', false);
    cmap = lines(nStimLevels);

    for chanInt = chanIntList

        % --- Plot 1: Mean EP at each stim level, pooled across all R blocks ---
        figure('Units','inches','Position',[1 1 8 5]);
        for condInt = 1:nStimLevels
            allTrials = [];
            for bi = 1:length(blocks)
                epochData = epochsEPblock{bi}{condInt};
                allTrials = cat(2, allTrials, squeeze(epochData(:, chanInt, :)));
            end
            meanEP = 1e6 * mean(allTrials, 2);
            plot(tEpoch, meanEP, 'linewidth', 2, 'color', cmap(condInt,:));
            hold on
        end
        xlim([-10 70])
        ylim([-350 350])
        ylabel('Voltage (\muV)')
        xlabel('time (ms)')
        legend(stimLabels, 'Location', 'best')
        title(sprintf('%s - Mean EP by stim level (all R blocks pooled), Chan %d', sid, chanInt))
        set(gca, 'fontsize', 12)
        line1 = vline(tBegin,'k:'); set(line1,'tag','vline','handlevisibility','off')
        line2 = vline(tEnd,'k:'); set(line2,'tag','vline','handlevisibility','off')

        % --- Plot 2: Same with CI shading ---
        figure('Units','inches','Position',[1 1 8 5]);
        for condInt = 1:nStimLevels
            allTrials = [];
            for bi = 1:length(blocks)
                epochData = epochsEPblock{bi}{condInt};
                allTrials = cat(2, allTrials, squeeze(epochData(:, chanInt, :)));
            end
            plotBTLError(tEpoch(:)', 1e6*allTrials, 'CI', cmap(condInt,:)');
        end
        xlim([-10 70])
        ylim([-350 350])
        ylabel('Voltage (\muV)')
        xlabel('time (ms)')
        h = flipud(findobj(gca,'Type','line'));
        legend(h, stimLabels, 'Location', 'best')
        title(sprintf('%s - Mean EP +/- 95%% CI by stim level, Chan %d', sid, chanInt))
        set(gca, 'fontsize', 12)
        line1 = vline(tBegin,'k:'); set(line1,'tag','vline','handlevisibility','off')
        line2 = vline(tEnd,'k:'); set(line2,'tag','vline','handlevisibility','off')

        % --- Plot 3: One figure per block, all stim levels overlaid ---
        for bi = 1:length(blocks)
            figure('Units','inches','Position',[1 1 8 5]);
            for condInt = 1:nStimLevels
                epochData = epochsEPblock{bi}{condInt};
                meanEP = 1e6 * mean(squeeze(epochData(:, chanInt, :)), 2);
                plot(tEpoch, meanEP, 'linewidth', 2, 'color', cmap(condInt,:));
                hold on
            end
            xlim([-10 70])
            ylim([-350 350])
            ylabel('Voltage (\muV)')
            xlabel('time (ms)')
            legend(stimLabels, 'Location', 'best')
            title(sprintf('%s - Block %d: %s, Chan %d', sid, blocks(bi), legendText{bi}, chanInt))
            set(gca, 'fontsize', 12)
            line1 = vline(tBegin,'k:'); set(line1,'tag','vline','handlevisibility','off')
            line2 = vline(tEnd,'k:'); set(line2,'tag','vline','handlevisibility','off')
        end
    end
end
