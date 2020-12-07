% Maps of individual cells� place fields along the track, for cells that were monitored together within a single imaging session.
fig_name = 'f1b_sortedRasters';

rng(0);

sm = SessManager;
d = sm.cons_usable(SessManager.special_sessions_usable_index('Mouse2022'));

right_mean_activity = mean( d.data_tensor(:,:,d.tr_dir== 1) , 3);
 left_mean_activity = mean( d.data_tensor(:,:,d.tr_dir==-1) , 3);

[~, right_peak_bin] = max(right_mean_activity, [], 2);
[~,  left_peak_bin] = max( left_mean_activity, [], 2);

[~, right_ordering] = sort(right_peak_bin);
[~,  left_ordering] = sort( left_peak_bin);

right_mean_activity = right_mean_activity(right_ordering, :);
 left_mean_activity =  left_mean_activity( left_ordering, :);

figure;
subplot(1,2,1);
imagesc(right_mean_activity, [0 1]);
xlabel '(running right)'
ylabel 'cell number'
ax = gca;
ax.XTick = [0.5 10.5 20.5];
ax.XTickLabel = {'0', 'L/2', 'L'};

hc = colorbar;
hc.Label.String = '\DeltaF/F_0 (Clamped to [0-1])';
hc.Label.Rotation = 270;
hc.Label.Position(1) = 4;
hc.Visible = false;

subplot(1,2,2);
imagesc(left_mean_activity, [0 1]);
xlabel '(running left)'
ax = gca;
ax.XTick = [0.5 10.5 20.5];
ax.XTickLabel = {'0', 'L/2', 'L'};

hc = colorbar;
hc.Label.String = '\DeltaF/F_0 (Clamped to [0-1])';
hc.Label.Rotation = 270;
hc.Label.Position(1) = 4;

suptitle(sprintf('simultaneously recorded cells\nordered by RFs'' peak-locat.'));

render_fig(fig_name);
% if exist(fig_name, 'dir') == 0, mkdir(fig_name); end
% savefig(fullfile(fig_name, [fig_name '.fig']));
% print(gcf, '-dpdf', fullfile(fig_name, [fig_name '.pdf']));