function multi_figure_format(varargin)
p = inputParser;
p.addOptional('fix_expo', {}, @iscell);
p.parse(varargin{:});

my_font = 'Helvetica LT Std';
my_font_size = 6*6/6.4;
my_line_width = 0.5;
my_boxsize = [0.8125 0.585].*0.98;

axs = get(gcf, 'Children');
for i = 1:numel(axs)
    ax = axs(i);
    ax.FontName = my_font;
    ax.FontSize = my_font_size;
    ax.Box = 'off';
    ax.LineWidth = my_line_width;
    ax.Units = 'inches';
    ax.Position = [ax.Position(1:2), my_boxsize];
    ax.TickLength = [0.02 0.02];
    ax.XLabel.FontName = my_font;
    ax.YLabel.FontName = my_font;
    ax.XLabel.FontSize = my_font_size;
    ax.YLabel.FontSize = my_font_size;
    if ~isempty(p.Results.fix_expo)
        Utils.fix_exponent(ax, p.Results.fix_expo{:});
    end
    for j = 1:numel(ax.Children)
        child = ax.Children(j);
        if strcmp(child.Type, 'text')
            child.FontName = my_font;
            child.FontSize = my_font_size;
        end
    end
end

end