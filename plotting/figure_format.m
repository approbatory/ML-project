function figure_format(varargin)
p = inputParser;
p.addOptional('boxsize', [0.8125 0.585].*0.98, @(x) isnumeric(x) && (numel(x) == 2));
p.addOptional('lwid', 0.5, @isscalar);
p.addOptional('fontsize', 6, @isscalar);
p.parse(varargin{:});

hf = gcf;
hf.Units = 'inches';
hf.Position = [hf.Position(1:2), p.Results.boxsize.*1.6];

ax = gca;
ax.FontSize = p.Results.fontsize;
ax.FontName = 'Arial';
ax.Box = 'off';
ax.LineWidth = p.Results.lwid;
ax.Units = 'inches';
ax.Position = [ax.Position(1:2), p.Results.boxsize];
ax.TickLength = [0.02 0.02];

for i = 1:numel(ax.Children)
    if strcmp(ax.Children(i).Type, 'text')
        ax.Children(i).FontSize = p.Results.fontsize;
    end
end
% box = [0 0 1.3 0.8];
%
% set(fig, 'DefaultLineLineWidth', 1, ...
%            'DefaultAxesFontName', 'Arial', ...
%            'DefaultAxesLineWidth', 1, ...
%            'DefaultAxesFontSize', 8, ...
%            'DefaultAxesBox', 'off', ...
%            'DefaultAxesColor', [1, 1, 1], ...
%            'DefaultFigureColor', [1, 1, 1], ...
%            'Units', 'inches',...
%            'DefaultFigureInvertHardcopy', 'off', ...
%            'DefaultFigurePaperUnits', 'inches', ...
%            'DefaultFigureUnits', 'inches', ...
%            'DefaultFigurePaperPosition', box, ...
%            'DefaultFigurePosition', box, ...
%            'DefaultFigurePaperSize', box(3:4),...
%            'Position', box,...
%            'Color', [1,1,1]);
% set(fig.CurrentAxes, 'Units', 'inches', 'Position', box/3);
%if before
%    set(gcf, 'Units', 'inches', 'Position', [0 0 1.3 1],...
%        'DefaultAxesFontSize', 8, 'DefaultAxesFontName', 'Arial');
%else
%    set(gca, 'Box', 'off');
%
%end
if false
    %%
    figure;
    fwid = 6;
    lwid = 0.5;
    plot((1:30)./30.*24);
    tx = xlabel('Time interval (d)');
    ty = ylabel('Median\newlineerror (cm)');
    xlim([-2.5 35]);
    ylim([0 24]);
    plt = Plot();
    plt.ShowBox = 'off';
    plt.BoxDim = [0.8125 0.585].*0.98;
    plt.FontName = 'Arial';
    plt.FontSize = fwid;
    %plt.XLabel = 'Time interval (d)';
    %plt.YLabel = 'Median\newlineerror (cm)';
    plt.LineWidth = lwid;
    plt.AxisLineWidth = lwid;
    plt.YMinorTick = 'off';
    plt.XMinorTick = 'off';
    plt.XTick = 0:5:35;
    plt.XTickLabel(2:2:end) = {[]};
    plt.YTick = 0:6:24;
    plt.YTickLabel(2:end-1) = {[]};
    text(3, 5, 'Example', 'FontSize', fwid, 'Color', 'red');
end
if false
    %%
    fwid = 6;
    lwid = 0.5;
    boxsize = [0.8125 0.585].*0.98;
    
    figure;
    plot((1:30)./30.*24);
    xlabel('Time interval (d)');
    ylabel('Median\newlineerror (cm)');
    xlim([-2.5 35]);
    ylim([0 24]);
    text(3, 5, 'Example', 'FontSize', fwid, 'Color', 'red');
    
    
    
    hf = gcf;
    hf.Units = 'inches';
    hf.Position = [hf.Position(1:2), boxsize.*1.6];
    
    ax = gca;
    ax.FontSize = fwid;
    ax.FontName = 'Arial';
    ax.Box = 'off';
    ax.LineWidth = lwid;
    ax.Units = 'inches';
    ax.Position = [ax.Position(1:2), boxsize];
    %ax.XTick = 0:5:35;
    %ax.XTickLabel(2:2:end) = {''};
    %ax.YTick = 0:6:24;
    %ax.YTickLabel(2:end-1) = {''};
    
end
end