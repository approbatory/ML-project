function ResizeFigure(figureHandle, newWidth, newHeight, units)
% ResizeFigure(figureHandle, newWidth, newHeight, units)
%
% Resizes the specified figure to have the specified width and
% height, specified in the given units.

set(figureHandle,'Units',units);
position = get(figureHandle, 'Position');
set(figureHandle, 'Position',[position(1:2) newWidth, newHeight]);



