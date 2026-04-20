% Reproduces the a23ed_chan_5_mean_measure.fig / confInt_measure.fig style
% for any list of (subject, channel, blocks, legendText) cases.
%
% Per case, at each stim level in stimInterest, emits:
%   {sid}_chan_{ch}_mean_measure_stim{s}.png   (and .fig)
%   {sid}_chan_{ch}_confInt_measure_stim{s}.png (and .fig)
%
% Each figure overlays the configured blocks (cbrewer Dark2) on a single
% panel and marks the baseline peak / trough times with heavy dashed
% vertical lines (vline + LW 2, like the reference .fig). Peak / trough
% come from pkLocsBlock / trLocsBlock, which extract_PP_peak_to_peak
% computes on the trial-averaged waveform, so exactly one peak line and
% one trough line are drawn per panel.
%
% Plotting is inline (not via analyze_EP_compare_multiple_blocks) to use
% the actual stim level being plotted for peak/trough lookup; that helper
% hardcodes condInt = 4.

close all; clear; clc
Z_ConstantsDBS_PairedPulse;
matlab_dir = 'MATLAB_Converted';
experiment = 'EP_Measurement';

% Analysis flags (match master_script_analyze_EP.m defaults)
avgTrials      = 0;
savePlot       = 0;
screenBadChans = 0;
plotCondAvg    = 0;
plotPkTr       = 0;
tryArtifact    = 0;
chanReref      = 1;
rerefMode      = 'none';
numAvg         = 3;

% --- USER CONFIG ------------------------------------------------------
% Stim levels (condInt indices into stimLevelUniq). Reference uses [4].
stimInterest = [3 4];

% Cases: sid | channel | blocks | legendText
% Defaults: a23ed chan 5 (reference), 41a73 chan 5 + 68574 chan 7 (matching).
cases = {
    'a23ed', 5, [2 3 6 8], {'baseline 2', 'post A/B 200 ms (5 minutes)', ...
                            'post A/B 200 ms (15 mins)', 'post A/A 200 ms (15 mins)'};
    'a23ed', 5, [5 6],     {'baseline 4 (matched pre for 15-min A/B 200)', 'post A/B 200 ms (15 mins)'};
    '41a73', 5, [5 6],     {'baseline 4 (pre A/A 200)', 'post A/A 200 ms'};
    '68574', 7, [6 7],     {'baseline 4 (pre A/A 200)', 'post A/A 200 ms'};
    '9f852', 4, [2 3],     {'baseline 2 (pre conditioning)', 'post A/B 25 ms'};
    };

% Plot limits and vline style (match a23ed_chan_5_mean_measure.fig)
xLimMs     = [-10 70];
yLimUV     = [-300 300];
xTicks     = [0 20 40 60];
yTicks     = [-200 -100 0 100 200];

% Per-subject overrides (smaller-amplitude subjects get tighter YLim/YTicks)
yLimPerSid = containers.Map;
yLimPerSid('68574') = [-140 140];
yTicksPerSid = containers.Map;
yTicksPerSid('68574') = [-100 0 100];
axFontSize = 22;
vlineLW    = 2;
figPosIn   = [1 1 7 7];  % square figure; axes forced square via pbaspect

saveDir = OUTPUT_DIR;
% --- END USER CONFIG --------------------------------------------------

if ~exist(saveDir, 'dir'); mkdir(saveDir); end

dataCache = containers.Map;

for c = 1:size(cases, 1)
    sid         = cases{c, 1};
    chanInt     = cases{c, 2};
    blocks      = cases{c, 3};
    legendText  = cases{c, 4};
    chanIntList = chanInt;  % consumed by prepare_EP_blocks

    assert(length(blocks) == length(legendText), ...
        'Case %d (%s): blocks and legendText must match in length', c, sid);

    cacheKey = sprintf('%s_%s', sid, strjoin(string(blocks), '_'));
    if ~isKey(dataCache, cacheKey)
        fprintf('\n=== Preparing %s (blocks: %s) ===\n', sid, num2str(blocks));
        prepare_EP_blocks;
        dataCache(cacheKey) = struct( ...
            'epochsEPblock', {epochsEPblock}, ...
            'tEpoch',        tEpoch, ...
            'stimLevelUniq', stimLevelUniq, ...
            'tBegin',        tBegin, ...
            'tEnd',          tEnd, ...
            'ECoGfs',        ECoGfs, ...
            'pkLocsBlock',   {pkLocsBlock}, ...
            'trLocsBlock',   {trLocsBlock});
    end
    d = dataCache(cacheKey);

    cmap = cbrewer('qual', 'Dark2', max(3, length(blocks)));
    cmap = cmap(1:length(blocks), :);
    cmap(cmap > 1) = 1; cmap(cmap < 0) = 0;

    tBeginSamp = d.ECoGfs * d.tBegin / 1e3;
    tEndSamp   = d.ECoGfs * d.tEnd   / 1e3;
    tPk        = (tBeginSamp:tEndSamp) / d.ECoGfs;

    if isKey(yLimPerSid, sid)
        currYLim   = yLimPerSid(sid);
        currYTicks = yTicksPerSid(sid);
    else
        currYLim   = yLimUV;
        currYTicks = yTicks;
    end

    for sIdx = 1:length(stimInterest)
        condInt = stimInterest(sIdx);
        if condInt > length(d.stimLevelUniq)
            fprintf('  %s: stim level %d unavailable (max %d); skipping\n', ...
                sid, condInt, length(d.stimLevelUniq));
            continue
        end
        stimUA = d.stimLevelUniq(condInt);

        % Peak / trough times from baseline block (first in blocks list),
        % computed on the trial-averaged waveform by extract_PP_peak_to_peak.
        pkIdx = d.pkLocsBlock{1}(chanInt, condInt);
        trIdx = d.trLocsBlock{1}(chanInt, condInt);
        pkMs  = 1e3 * tPk(pkIdx);
        trMs  = 1e3 * tPk(trIdx);

        % --- Mean figure ---
        meanFig = figure('Units', 'inches', 'Position', figPosIn);
        hold on
        hBlocks = gobjects(1, length(blocks));
        for i = 1:length(blocks)
            trials = squeeze(d.epochsEPblock{i}{condInt}(:, chanInt, :));
            hBlocks(i) = plot(d.tEpoch, 1e6 * mean(trials, 2), ...
                'LineWidth', 2, 'Color', cmap(i, :));
        end
        xlim(xLimMs); ylim(currYLim);
        xlabel('time (ms)'); ylabel('Voltage (\muV)');
        title(sprintf('Channel %d, %d \\muA (mean)', chanInt, stimUA));
        set(gca, 'FontSize', axFontSize, 'XTick', xTicks, 'YTick', currYTicks, ...
            'TitleFontSizeMultiplier', 0.8, 'LabelFontSizeMultiplier', 1.0);
        pbaspect([1 1 1]);  % force axes box to be square
        legend(hBlocks, legendText);
        pkLine = vline(pkMs, 'k--'); set(pkLine, 'LineWidth', vlineLW, 'Tag', 'vline');
        trLine = vline(trMs, 'k--'); set(trLine, 'LineWidth', vlineLW, 'Tag', 'vline');

        blockStr = strjoin(string(blocks), '_');
        meanName = sprintf('%s_chan_%d_blocks_%s_mean_measure_stim%d', sid, chanInt, blockStr, condInt);
        exportgraphics(meanFig, fullfile(saveDir, [meanName '.png']), 'Resolution', 600);
        exportgraphics(meanFig, fullfile(saveDir, [meanName '.eps']), 'ContentType', 'vector');
        savefig(meanFig, fullfile(saveDir, [meanName '.fig']));

        % --- Confidence interval figure ---
        confFig = figure('Units', 'inches', 'Position', figPosIn);
        hold on
        for i = 1:length(blocks)
            trials = squeeze(d.epochsEPblock{i}{condInt}(:, chanInt, :));
            plotBTLError(d.tEpoch(:)', 1e6 * trials, 'CI', cmap(i, :)');
        end
        xlim(xLimMs); ylim(currYLim);
        xlabel('time (ms)'); ylabel('Voltage (\muV)');
        title(sprintf('Channel %d, %d \\muA (95%% CI)', chanInt, stimUA));
        set(gca, 'FontSize', axFontSize, 'XTick', xTicks, 'YTick', currYTicks, ...
            'TitleFontSizeMultiplier', 0.8, 'LabelFontSizeMultiplier', 1.0);
        pbaspect([1 1 1]);  % force axes box to be square
        hAll = flipud(findobj(gca, 'Type', 'line'));
        if length(hAll) >= length(blocks)
            legend(hAll(1:length(blocks)), legendText);
        end

        confName = sprintf('%s_chan_%d_blocks_%s_confInt_measure_stim%d', sid, chanInt, blockStr, condInt);
        exportgraphics(confFig, fullfile(saveDir, [confName '.png']), 'Resolution', 600);
        exportgraphics(confFig, fullfile(saveDir, [confName '.eps']), 'ContentType', 'vector');
        savefig(confFig, fullfile(saveDir, [confName '.fig']));

        fprintf('  %s chan %d stim %d (%.1f mA): pk=%.2f ms tr=%.2f ms  [mean + confInt saved]\n', ...
            sid, chanInt, condInt, stimUA/1e3, pkMs, trMs);
    end
end

fprintf('\nDone. %d case(s) plotted.\n', size(cases, 1));
