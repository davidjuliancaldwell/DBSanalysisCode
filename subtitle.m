function [ax,h]=subtitle(text)
%
%Centers a title over a group of subplots.
%Returns a handle to the title and the handle to an axis.
% [ax,h]=subtitle(text)
%           returns handles to both the axis and the title.
% ax=subtitle(text)
%           returns a handle to the axis only.

% increase the first two values in the axes label to further offset the
% title
ax=axes('Units','Normal','Position',[.1 .1 .85 .85],'Visible','off');
set(get(ax,'Title'),'Visible','on')
title(text, 'FontSize', 20); % modified DJC 9-8-2015 for bigger size 
if (nargout < 2)
    return
end
h=get(ax,'Title');

end