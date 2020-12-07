function render_fig(fig_name)
if exist(fig_name, 'dir') == 0, mkdir(fig_name); end
savefig(fullfile(fig_name, [fig_name '.fig']));
print(gcf, '-dpdf', fullfile(fig_name, [fig_name '.pdf']));