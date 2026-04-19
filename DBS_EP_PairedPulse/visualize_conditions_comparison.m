% Standalone driver for "comparison of conditions, average plot" figures.
% Overlays all configured blocks at each stim level (3 and 4 by default)
% for the given channel(s), producing mean-only and mean + 95% CI panels.
%
% Wraps the existing prepare_EP_blocks + analyze_EP_compare_multiple_blocks
% pipeline. Edit the USER CONFIG block for any (subject, channel, blocks)
% combination; the defaults plot 9f852 channel 7 across all its paired
% pulse blocks.
%
% Saves figures via exportgraphics (SaveFig has a macOS path-prefix bug
% that writes to c:/Tim/... on Unix systems).

close all; clear; clc
Z_ConstantsDBS_PairedPulse;
matlab_dir = 'MATLAB_Converted';
experiment = 'EP_Measurement';

% Analysis flags (match master_script_analyze_EP.m defaults)
avgTrials      = 0;
savePlot       = 0;   % keep 0 -- exportgraphics saves below
screenBadChans = 0;
plotCondAvg    = 0;
plotPkTr       = 0;
tryArtifact    = 0;
chanReref      = 1;
rerefMode      = 'none';
numAvg         = 3;

% --- USER CONFIG ------------------------------------------------------
% Subject ID and channel(s) to plot
sid         = '9f852';
chanIntList = [7];

% Blocks to overlay and their legend labels (must be same length).
% Default values taken from master_script_analyze_EP.m for 9f852.
blocks = [2 3 4 5 6 7 10 11 12];
legendText = {'baseline 2 (pre conditioning)', 'post A/B 25 ms', ...
    'baseline 3 (post 25 ms)', 'post A/A 25 ms', ...
    'baseline 4 (post 25 ms A/A)', 'post A/B 200 ms', ...
    'baseline 5 - post A/B 200 ms 12 min later', ...
    'post A/B 25 ms second time', 'baseline 6'};

% Where to write PNGs; defaults to the OUTPUT_DIR from Z_Constants.
saveDir = OUTPUT_DIR;
% --- END USER CONFIG --------------------------------------------------

assert(length(blocks) == length(legendText), ...
    'blocks and legendText must be same length');

fprintf('\n=== %s, chan %s, %d blocks ===\n', sid, mat2str(chanIntList), length(blocks));
fprintf('Preparing EP blocks...\n');
prepare_EP_blocks

fprintf('Running analyze_EP_compare_multiple_blocks...\n');
figsBefore = findall(groot, 'Type', 'figure');
analyze_EP_compare_multiple_blocks
figsAfter  = findall(groot, 'Type', 'figure');

% --- Save new figures via exportgraphics ---
newFigs = setdiff(figsAfter, figsBefore);
[~, ord] = sort(arrayfun(@(f) f.Number, newFigs));
newFigs = newFigs(ord);

if ~exist(saveDir, 'dir'); mkdir(saveDir); end
chanStr = strjoin(arrayfun(@num2str, chanIntList, 'UniformOutput', false), '_');

fprintf('\nSaving %d figure(s) to %s\n', length(newFigs), saveDir);
for k = 1:length(newFigs)
    f  = newFigs(k);
    ax = findobj(f, 'Type', 'axes');
    ttlstr = '';
    if ~isempty(ax)
        t = get(ax(1), 'Title');
        if ~isempty(t) && ~isempty(t.String)
            rawTitle = t.String;
            if iscell(rawTitle); rawTitle = strjoin(rawTitle, ' '); end
            ttlstr = regexprep(rawTitle, '[^a-zA-Z0-9]+', '_');
            ttlstr = regexprep(ttlstr, '^_+|_+$', '');
        end
    end
    fname = sprintf('conditions_%s_chan%s_fig%02d_%s.png', ...
        sid, chanStr, k, ttlstr);
    exportgraphics(f, fullfile(saveDir, fname), 'Resolution', 600);
    fprintf('  [%d/%d] %s\n', k, length(newFigs), fname);
end

fprintf('\nDone.\n');
