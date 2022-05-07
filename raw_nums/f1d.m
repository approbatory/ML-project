function f1d
%% Figure 1d
% Raw numbers for:
% - LineX
% - RightRealY, RightRealShadedY
% - LeftRealY, LeftRealShadedY
% - RightShuffledY, RightShuffledShadedY
% - LeftShuffledY, LeftShuffledShadedY
rng(0);

org = Org;
org.load_definitions;

data = plot_noise_corr_by_dist(org);
make_xlsx(data, 'f1d');
end

%% helper functions
function data = plot_noise_corr_by_dist(org)
figure;
ncd_r = cell2mat(org.fetch('nc_by_dist_r')');
r_h = serrorbar(0:19, mean(ncd_r), sem(ncd_r).*1.96, 'lineprops', {'-', 'Color', [0 0.5 0]});
hold on;

data.LineX = 0:19;
data.RightRealY = mean(ncd_r);
data.RightRealShadedY = sem(ncd_r).*1.96;

ncd_l = cell2mat(org.fetch('nc_by_dist_l')');
l_h = serrorbar(0:19, mean(ncd_l), sem(ncd_l).*1.96, 'lineprops', {'-', 'Color', [0.8 0.5 0]});

data.LeftRealY = mean(ncd_l);
data.LeftRealShadedY = sem(ncd_l).*1.96;

ncd_r = cell2mat(org.fetch('nc_by_dist_r_shuf')');
r_hs = serrorbar(0:19, mean(ncd_r), sem(ncd_r).*1.96, 'lineprops', {':', 'Color', 1/4*[1 1 1]});% [0 0.5 0] * 0.75
hold on;

data.RightShuffledY = mean(ncd_r);
data.RightShuffledShadedY = sem(ncd_r).*1.96;

ncd_l = cell2mat(org.fetch('nc_by_dist_l_shuf')');
l_hs = serrorbar(0:19, mean(ncd_l), sem(ncd_l).*1.96, 'lineprops', {':', 'Color', 1/4*[1 1 1]});% [0.8 0.5 0] * 0.75

data.LeftShuffledY = mean(ncd_l);
data.LeftShuffledShadedY = sem(ncd_l).*1.96;

xlabel 'RFs'' peak distance (cm)'
ylabel 'Pairwise noise correlations'
legend([r_h.mainLine l_h.mainLine], 'Running right', 'Running left');
legend boxoff

ylim([-0.01 inf]);
xlim([0 20]);

text(20, -0.003, 'Shuffled', 'Color', [1 1 1]/4, 'HorizontalAlignment', 'right');

figure_format([1.3 1.3]);
%set(gca, 'XAxisLocation', 'origin');
%set(gca, 'YAxisLocation', 'origin');
set(gca, 'XTick', [0 10 20]);
set(gca, 'XTickLabel', {'0', '60', '120'});
end