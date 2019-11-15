function h = HorizontalLine(ax, y, varargin)
% h = HorizontalLine(ax, y, varargin)
% 
% Inserts a horizontal line with at the specified x value in the
% specified axis. The line is red by default, but optional arguments
% can set its style.
%
% H is a handle to the line.

hold(ax,'on');
xl = xlim(ax);
if (isempty(varargin))
  h = plot(xl,[y y],'r');
else
  h = plot(xl,[y,y],varargin{:});
end
