function figure_format(varargin)
p = inputParser;
p.addOptional('boxsize', [0.8125 0.585].*0.98, @(x) isnumeric(x) && (numel(x) == 2));
p.addParameter('factor', 1, @isscalar);
p.addParameter('fig_factor', 1.6, @isscalar);
p.addParameter('lwid', 0.5, @isscalar);
p.addParameter('fontsize', 6, @isscalar);
p.parse(varargin{:});

my_font = 'Helvetica';

hf = gcf;
hf.Units = 'inches';
hf.Position = [hf.Position(1:2), p.Results.boxsize.*p.Results.fig_factor.*p.Results.factor];



ax = gca;
ax.FontName = my_font;
ax.FontSize = p.Results.fontsize;
ax.Box = 'off';
ax.LineWidth = p.Results.lwid;
ax.Units = 'inches';
ax.Position = [ax.Position(1:2), p.Results.boxsize.*p.Results.factor];
ax.TickLength = [0.02 0.02];
ax.YColor = 'k';
ax.XColor = 'k';

ax.XLabel.FontName = my_font;
ax.YLabel.FontName = my_font;
ax.XLabel.FontSize = p.Results.fontsize;
ax.YLabel.FontSize = p.Results.fontsize;
for i = 1:numel(ax.Children)
    if strcmp(ax.Children(i).Type, 'text')
        ax.Children(i).FontName = my_font;
        ax.Children(i).FontSize = p.Results.fontsize;
    end
end

ax_children = hf.Children;
for i = 1:numel(ax_children)
    if ~isequal(ax_children(i), ax)
        ax_children(i).Color = 'k';
    end
end