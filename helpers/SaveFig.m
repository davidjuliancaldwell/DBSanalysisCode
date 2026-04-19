function SaveFig(path,filename,filetype,opt)
% path     - directory to save into (absolute on any OS, or relative -> legacy fallback)
% filename - file name (without extension)
% filetype - 'png' (default), 'eps', 'jpg', or 'svg'
% opt      - optional resolution string (e.g. '-r600')
%
% Cross-platform absolute-path handling:
%   Windows: '<letter>:' prefix, e.g. 'C:\...'
%   Unix:    '/' prefix, e.g. '/Users/...'
% Relative paths fall through to a legacy default directory.
%
% Patched from MATLAB_ECoG_code/Visualization/SaveFig.m (original only
% recognized 'C:/' and 'D:/'; any Unix path was prefixed with
% 'c:/Tim/research/script/generated_figs/', producing literal 'c:' dirs
% on macOS).

if ~exist('filetype','var')
    filetype = 'png';
end

% Detect absolute paths across platforms
isAbsoluteWin = length(path) >= 2 && path(2) == ':' && isletter(path(1));
isAbsoluteUnix = ~isempty(path) && path(1) == '/';

if ~(isAbsoluteWin || isAbsoluteUnix)
    % Legacy fallback for relative paths (kept for backward compatibility)
    if path(end) ~= '/' && path(end) ~= '\'
        path(end+1) = '/';
    end
    path = ['c:/Tim/research/script/generated_figs/' path];
end

% Ensure trailing separator so TouchDir + filename concat works
if path(end) ~= '/' && path(end) ~= '\'
    path(end+1) = '/';
end
TouchDir(path);
destFile = [path filename '.'];

warning('off', 'MATLAB:print:adobecset:DeprecatedOption')
set(gcf,'PaperPositionMode','auto');

if exist('opt','var') && ~isempty(opt)
    switch filetype
        case 'jpg'
            print('-djpeg', '-noui','-cmyk', '-painters', [destFile filetype], opt);
        case 'png'
            print('-dpng',  '-noui', '-opengl',           [destFile filetype], opt);
        case 'eps'
            print('-dpsc2', '-noui', '-painters',         [destFile filetype], opt);
        case 'svg'
            print('-dsvg',  '-noui', '-painters',         [destFile filetype], opt);
    end
else
    switch filetype
        case 'jpg'
            print('-djpeg', '-noui','-cmyk', '-painters', [destFile filetype]);
        case 'png'
            print('-dpng',  '-noui', '-opengl',           [destFile filetype]);
        case 'eps'
            print('-dpsc2', '-noui', '-painters',         [destFile filetype]);
        case 'svg'
            print('-dsvg',  '-noui', '-painters',         [destFile filetype]);
    end
end
