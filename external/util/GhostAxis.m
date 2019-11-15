function GhostAxis(varargin)
% GhostAxis(varargin)
%
% Converts the axis AX (defaults to current axis) into a 'ghost axis'
% that is invisible but still appears in the figure when it's exported
% by turning the axis off, setting the color to GHOSTCOLOR (defaults
% to white), and plotting a single point in the middle with the color
% GHOSTCOLOR. Optional arguments 'xlabel' and 'ylabel' can be provided as
% key-value pairs with the corresponding arguments to be passed on to
% XLABEL and YLABEL.
%
% This function is useful when you want to make a multi-stage figure
% where the axis appear in turn. If you just the axis to invisbile, it
% won't be exported properly, so we just 'ghost' the axes. 
%
% Optional arguments:
%
% axis: The axis to ghost. Defaults to the current axis.
% ghostColor: The color to make the axis. Defaults to white.
% xlabel: Cell array containing input arguments for xlabel.
% ylabel: Cell array containing input arguments for ylabel.

p = inputParser;
p.addOptional('axis', gca);
p.addOptional('ghostColor',[1 1 1]);
p.addOptional('xlabel', []);
p.addOptional('ylabel', []);
p.parse(varargin{:});
opts = p.Results;

ax         = opts.axis;
ghostColor = opts.ghostColor;

if (~iscell(opts.xlabel))
  opts.xlabel = {opts.xlabel};
end

if (~iscell(opts.ylabel))
  opts.ylabel = {opts.ylabel};
end

xl = xlim(ax);
yl = ylim(ax);
plot(ax, xl(1)+diff(xl)/2, yl(1) + diff(yl)/2, 'Color',ghostColor);
set(ax,'Color',ghostColor);
set(ax,'xtick',[],'ytick',[],'xcolor',ghostColor,'ycolor',ghostColor);
h = xlabel(ax, opts.xlabel{:});
set(h,'Color', ghostColor);
h = ylabel(ax, opts.ylabel{:});
set(h, 'Color', ghostColor);
