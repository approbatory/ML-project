% Pairwise correlated fluctuations between neurons as a function of the distance between the RFs’ peaks (3,941,140 pairs from 110 sessions from 12 mice). There is a positive mean correlation between RFs with peaks separated by up to an ? of the track length.
fig_name = 'f1d_noiseCorrByPeakDistance';

org = Org;
org.load_definitions;
%%
figure;
ncd_r = cell2mat(org.fetch('nc_by_dist_r')');
r_h = serrorbar(0:19, mean(ncd_r), sem(ncd_r).*1.96, 'lineprops', {'-', 'Color', [0 0.5 0]});
hold on;
ncd_l = cell2mat(org.fetch('nc_by_dist_l')');
l_h = serrorbar(0:19, mean(ncd_l), sem(ncd_l).*1.96, 'lineprops', {'-', 'Color', [0.8 0.5 0]});

l_ = refline(0,0); l_.Color = 'k';
xlabel 'RFs'' peak distance'
ylabel 'Pairwise noise correlations'
legend([r_h.mainLine l_h.mainLine], 'Running right', 'Running left');
legend boxoff

figure_format([1.3 1.3]);
set(gca, 'XTick', [0 10 20]);
set(gca, 'XTickLabel', {'0', 'L/2', 'L'});

render_fig(fig_name);