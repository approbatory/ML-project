% to run all:
% S = dir('full_records'); for fs = S', if ~fs.isdir, o = Analyzer.recreate(fullfile(fs.folder, fs.name)); pls_origin; end, end; close all; clear; disp('done');
printit = true;

if printit
    mouse_id = split(o.res.source,'/');
    mouse_id = mouse_id{end};
    mouse_id = split(mouse_id, '-');
    mouse_id = [mouse_id{1} mouse_id{2}];
    
    fig_dir = ['graphs2/analyzer_figs_origin/large/' mouse_id];
    if ~exist(fig_dir, 'dir')
        mkdir(fig_dir);
    end
end
%% f_PLS_vis

[XL, ~, act,~,~,~,~,stats] = plsregress(o.data.X.fast, zscore([o.data.y.scaled, o.data.y.direction]), 2);
X_mean = mean(o.data.X.fast);

act_origin = -X_mean * stats.W;

ret = o.bin_data(false, false);
ret_shuf = o.bin_data(true, false);

f_PLS_vis = figure;
hold on;
scatter(act(:,1), act(:,2), 20, o.data.y.scaled,...
    'filled', 'MarkerFaceAlpha', 0.05);
%scatter(act_origin(1), act_origin(2), 20, 'r', 'filled');
xlabel PLS1; ylabel PLS2; title '2D PLS projections'


q_ = ret.bin_X.mean*stats.W;
for b = 1:size(q_,1)
    m_vec = plot([act_origin(1) q_(b,1)],...
        [act_origin(2) q_(b,2)], '-or');
end

for b = 1:size(q_,1)-2
    t_vec = plot([q_(b,1) q_(b+2,1)],...
        [q_(b,2) q_(b+2,2)], '-g');
end

fw_vars = ret.fw.latent(:,1);
bw_vars = ret.bw.latent(:,1);
q_f = q_(1:2:end,:);
z_f = fw_vars.*ret.fw.princ*stats.W + q_f;
mz_f = -fw_vars.*ret.fw.princ*stats.W + q_f;

q_b = q_(2:2:end,:);
z_b = bw_vars.*ret.bw.princ*stats.W + q_b;
mz_b = -bw_vars.*ret.bw.princ*stats.W + q_b;
for b = 1:o.opt.n_bins
    n_vec = plot([q_f(b,1) z_f(b,1)],...
        [q_f(b,2) z_f(b,2)], '-k');
    plot([q_b(b,1) z_b(b,1)],...
        [q_b(b,2) z_b(b,2)], '-k');
    plot([q_f(b,1) mz_f(b,1)],...
        [q_f(b,2) mz_f(b,2)], '-k');
    plot([q_b(b,1) mz_b(b,1)],...
        [q_b(b,2) mz_b(b,2)], '-k');
end
legend([m_vec, t_vec, n_vec],...
    'Mean activity', 'Curve tangent', 'Principal noise');
legend Location best
%% f_angles_by_bin
f_angles_by_bin = figure;
subplot(2,1,1); hold on;
plot(ret.fw.angles.noise_tangent, 'DisplayName', 'Noise/Tangent');
plot(ret.fw.angles.noise_mean, 'DisplayName', 'Noise/Mean');
plot(ret.fw.angles.mean_tangent, 'DisplayName', 'Mean/Tangent');
line(xlim,[90 90], 'Color', 'k', 'DisplayName', 'orthogonal');
l = legend('boxoff'); l.NumColumns = 2; legend Location best
title 'Dot product angles (forward direction)'
xlabel 'Bin index'
ylabel 'Angle (degrees)'

subplot(2,1,2); hold on; title 'bw'
plot(ret.bw.angles.noise_tangent, 'DisplayName', 'Noise/Tangent');
plot(ret.bw.angles.noise_mean, 'DisplayName', 'Noise/Mean');
plot(ret.bw.angles.mean_tangent, 'DisplayName', 'Mean/Tangent');
line(xlim,[90 90], 'Color', 'k', 'DisplayName', 'orthogonal');
l = legend('boxoff'); l.NumColumns = 2; legend Location best
title 'Dot product angles (backward direction)'
xlabel 'Bin index'
ylabel 'Angle (degrees)'
%%
h2_1 = figure;
subplot(2,1,1); hold on;
plot(ret_shuf.fw.angles.noise_tangent, 'DisplayName', 'Noise/Tangent');
plot(ret_shuf.fw.angles.noise_mean, 'DisplayName', 'Noise/Mean');
plot(ret_shuf.fw.angles.mean_tangent, 'DisplayName', 'Mean/Tangent');
line(xlim,[90 90], 'Color', 'k', 'DisplayName', 'orthogonal');
l = legend('boxoff'); l.NumColumns = 2; legend Location best
title 'Dot product angles, on shuffle (forward direction)'
xlabel 'Bin index'
ylabel 'Angle (degrees)'

subplot(2,1,2); hold on; title 'bw'
plot(ret_shuf.bw.angles.noise_tangent, 'DisplayName', 'Noise/Tangent');
plot(ret_shuf.bw.angles.noise_mean, 'DisplayName', 'Noise/Mean');
plot(ret_shuf.bw.angles.mean_tangent, 'DisplayName', 'Mean/Tangent');
line(xlim,[90 90], 'Color', 'k', 'DisplayName', 'orthogonal');
l = legend('boxoff'); l.NumColumns = 2; legend Location best
title 'Dot product angles, on shuffle (backward direction)'
xlabel 'Bin index'
ylabel 'Angle (degrees)'
%% f_angles_hist
f_angles_hist = figure;
subplot(3,1,1); hold on;
w_ = [ret.fw.angles.noise_tangent, ret.bw.angles.noise_tangent];
mw_ = mean(w_);
w_s = [ret_shuf.fw.angles.noise_tangent, ret_shuf.bw.angles.noise_tangent];
mw_s = mean(w_s);
histogram(w_s, 'DisplayName', 'Noise (shuffled)/Tangent', 'FaceColor', [0.8 0.8 0.9]);
histogram(w_, 'DisplayName', 'Noise/Tangent', 'FaceColor', [0 0 1]);
legend boxoff; legend Location west
line([mw_s mw_s], ylim, 'Color', 'k', 'LineStyle', '--', 'DisplayName', 'Mean angle (shuffled)');
line([mw_ mw_], ylim, 'Color', 'k', 'LineStyle', '-', 'LineWidth', 2, 'DisplayName', 'Mean angle');
xlim([0 90]);
xlabel 'Angle (degrees)'
ylabel 'Frequency'

subplot(3,1,2); hold on;
w_ = [ret.fw.angles.noise_mean, ret.bw.angles.noise_mean];
mw_ = mean(w_);
w_s = [ret_shuf.fw.angles.noise_mean, ret_shuf.bw.angles.noise_mean];
mw_s = mean(w_s);
histogram(w_s, 'DisplayName', 'Noise (shuffled)/Mean', 'FaceColor', [0.9 0.8 0.8]);
histogram(w_, 'DisplayName', 'Noise/Mean', 'FaceColor', [1 0 0]);
legend boxoff; legend Location west
line([mw_s mw_s], ylim, 'Color', 'k', 'LineStyle', '--', 'DisplayName', 'Mean angle (shuffled)');
line([mw_ mw_], ylim, 'Color', 'k', 'LineStyle', '-', 'LineWidth', 2, 'DisplayName', 'Mean angle');
xlim([0 90]);
xlabel 'Angle (degrees)'
ylabel 'Frequency'

subplot(3,1,3);
w_ = [ret.fw.angles.mean_tangent, ret.bw.angles.mean_tangent];
mw_ = mean(w_);
histogram(w_, 'DisplayName', 'Mean/Tangent', 'FaceColor', 'y');
legend boxoff; legend Location west
%line([90 90], ylim, 'Color', 'k', 'LineWidth', 2, 'DisplayName', 'Orthogonal');
line([mw_ mw_], ylim, 'Color', 'k', 'LineWidth', 2, 'DisplayName', 'Mean angle');
xlim([45 135]);
xlabel 'Angle (degrees)'
ylabel 'Frequency'
%% f_noise_spect
f_noise_spect = figure; hold on;
tot_lat = [ret.fw.latent; ret.bw.latent];
m_lat = mean(tot_lat);
s_lat = std(tot_lat);

tot_lat_shuf = [ret_shuf.fw.latent; ret_shuf.bw.latent];
m_lat_shuf = mean(tot_lat_shuf);
s_lat_shuf = std(tot_lat_shuf);
errorbar(m_lat, s_lat, 'DisplayName', 'unshuffled');
errorbar(m_lat_shuf, s_lat_shuf, 'DisplayName', 'shuffled');
set(gca, 'YScale', 'log');
legend
xlabel 'Principal component index'
ylabel 'Variance'
title 'Effect of shuffling on noise spectrum'
%% printing section
if printit
    printer = @(name) print('-dpng', fullfile(fig_dir, name));
    figure(f_PLS_vis); printer('PLS_vis.png');
    figure(f_angles_by_bin); printer('angles_by_bin.png');
    figure(f_angles_hist); printer('angles_hist.png');
    figure(f_noise_spect); printer('noise_spect.png');
end