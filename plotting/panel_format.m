function panel_format(ax)

my_font = 'Arial';
my_font_size = 6;
my_line_width = 0.5;

if ~exist('ax', 'var')
    ax = gca;
end
ax.FontName = my_font;
ax.FontSize = my_font_size;
ax.Box = 'off';
ax.LineWidth = my_line_width;
ax.TickLength = [0.02 0.02];
ax.YColor = 'k';
ax.XColor = 'k';

ax.XLabel.FontName = my_font;
ax.YLabel.FontName = my_font;
ax.XLabel.FontSize = my_font_size;
ax.YLabel.FontSize = my_font_size;
ax.Title.FontName = my_font;
ax.Title.FontSize = my_font_size;

for i = 1:numel(ax.Children)
    if strcmp(ax.Children(i).Type, 'text')
        ax.Children(i).FontName = my_font;
        ax.Children(i).FontSize = my_font_size;
    end
end
