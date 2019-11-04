function h = VerticalLine(ax, x, varargin)
% h = VerticalLine(ax, x, varargin)
% 
% Inserts a vertical line with at the specified x value in the
% specified axis. The line is rd by default, but optional arguments
% can set its style. A handle to the line is returned.

hold(ax,'on');
yl = ylim(ax);
if (isempty(varargin))
  h = plot([x x],yl,'r');
else
  h = plot([x x],yl,varargin{:});
end
