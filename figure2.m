%creating figure 2
svg_save_dir = 'figure2_svg';
print_svg = @(name) print('-dsvg', fullfile(svg_save_dir, [name '.svg']));
print_png = @(name) print('-dpng', '-r1800', fullfile(svg_save_dir, [name '.png']));
%% PLS representations
m_i = 4;
use_zscore = true;
[source_path, mouse_name] = DecodeTensor.default_datasets(m_i);
opt = DecodeTensor.default_opt;
opt.restrict_trials = -1;

[T, d] = DecodeTensor.tensor_loader(source_path, mouse_name, opt);
% opt.neural_data_type = 'FST_events';
% [T_fst, d_fst] = DecodeTensor.tensor_loader(source_path, mouse_name, opt);
% opt.neural_data_type = 'spikeDeconv';
% [T_oasis, d_oasis] = DecodeTensor.tensor_loader(source_path, mouse_name, opt);

[X, ks] = DecodeTensor.tensor2dataset(T, d);

% [X_fst, ks_fst] = DecodeTensor.tensor2dataset(T_fst, d_fst);
% [X_oasis, ks_oasis] = DecodeTensor.tensor2dataset(T_oasis, d_oasis);
X_shuf = shuffle(X, ks);
if use_zscore
    X = zscore(X);
    X_shuf = zscore(X_shuf);
end
% X_fst_shuf = shuffle(X_fst, ks_fst);
% X_oasis_shuf = shuffle(X_oasis, ks_oasis);

[XS, stats, origin] = Utils.pls_short(X, [ceil(ks/2), mod(ks,2)]);
figure;
scatter(XS(:,1), XS(:,2), 4, 'b','filled', 'MarkerFaceAlpha', 0.02); hold on;
scatter(origin(1), origin(2), 10, 'k');
axis equal;
xlabel PLS1; ylabel PLS2;
%xlim(xlim.*0.8);
%ylim(ylim-0.01);
xl_ = xlim;
yl_ = ylim;
set(gca, 'XTick', []);
set(gca, 'YTick', []);
figure_format;
print_png('PLS_unshuf');

figure;
XS_s = (X_shuf - mean(X_shuf)) * stats.W;
origin_s = - mean(X_shuf) * stats.W;
scatter(XS_s(:,1), XS_s(:,2), 4, 'r','filled', 'MarkerFaceAlpha', 0.02); hold on;
scatter(origin_s(1), origin_s(2), 10, 'k');
xlabel PLS1; ylabel PLS2;
axis equal;
xlim(xl_); ylim(yl_);
set(gca, 'XTick', []);
set(gca, 'YTick', []);
figure_format;
print_png('PLS_shuf');

% [~, stats] = Utils.pls_plot(X_fst, [ceil(ks_fst/2), mod(ks_fst,2)]);
% suptitle('Using FST\_events, unshuffled');
% colormap jet;
% xl_ = xlim;
% yl_ = ylim;
% 
% Utils.pls_plot(X_fst_shuf, [ceil(ks_fst/2), mod(ks_fst,2)], stats, xl_, yl_);
% suptitle('Using FST\_events, shuffled');
% colormap jet;
% 
% [~, stats] = Utils.pls_plot(X_oasis, [ceil(ks_oasis/2), mod(ks_oasis,2)]);
% suptitle('Using spikeDeconv, unshuffled');
% colormap jet;
% xl_ = xlim;
% yl_ = ylim;
% 
% Utils.pls_plot(X_oasis_shuf, [ceil(ks_oasis/2), mod(ks_oasis,2)], stats, xl_, yl_);
% suptitle('Using spikeDeconv, shuffled');
% colormap jet;
% %% Measure noise in direction of signal % this section has been generalized in DecodeTensor.noise_properties
% %bin 10 as an example
% 
% X10 = squeeze(T(:,10,d==1)).';
% X_noise = X10 - mean(X10);
% Noise_Cov = cov(X_noise);
% 
% X11 = squeeze(T(:,11,d==1)).';
% X9 = squeeze(T(:,9,d==1)).';
% Signal_Direction_pre = mean(X10) - mean(X9);
% Signal_Direction_post = mean(X11) - mean(X10);
% Signal_Direction = (Signal_Direction_post + Signal_Direction_pre)/2;
% Signal_Direction = Signal_Direction.' ./ norm(Signal_Direction);
% 
% Noise_in_Direction_of_Signal = Signal_Direction.' * Noise_Cov * Signal_Direction;
% [~, latent] = pcacov(Noise_Cov);
% 
% 
% T_s = DecodeTensor.shuffle_tensor(T, d);
% X10 = squeeze(T_s(:,10,d==1)).';
% X_noise = X10 - mean(X10);
% Noise_Cov = cov(X_noise);
% 
% X11 = squeeze(T_s(:,11,d==1)).';
% X9 = squeeze(T_s(:,9,d==1)).';
% Signal_Direction_pre = mean(X10) - mean(X9);
% Signal_Direction_post = mean(X11) - mean(X10);
% Signal_Direction = (Signal_Direction_post + Signal_Direction_pre)/2;
% Signal_Direction = Signal_Direction.' ./ norm(Signal_Direction);
% 
% Noise_in_Direction_of_Signal_s = Signal_Direction.' * Noise_Cov * Signal_Direction;
% [~, latent_s] = pcacov(Noise_Cov);
% %% Noise spectrum changes, noise loadings, and total noise in signal direction
% res = DecodeTensor.noise_properties(T, d, false);
% figure; 
% subplot(1,2,1);
% plot(mean(cell2mat(res.el_pre(:)).^2));
% hold on;
% plot(mean(cell2mat(res.el_rnd(:)).^2));
% ylim([0 0.06]);
% subplot(1,2,2);
% plot(mean(cell2mat(res.el_pre_s(:)).^2));
% hold on;
% plot(mean(cell2mat(res.el_rnd_s(:)).^2));
% ylim([0 0.06]);
% 
% 
% figure;
% subplot(1,2,1);
% loglog(mean(cell2mat(res.noise_spectrum(:).'),2));
% ylim([1e-3 1e3]);
% subplot(1,2,2);
% loglog(mean(cell2mat(res.noise_spectrum_s(:).'),2));
% ylim([1e-3 1e3]);
% 
% res_corr = DecodeTensor.noise_properties(T, d, true);
% 
% figure;
% sh_by_unsh = cell2mat(res.nd_pre(:)) ./ cell2mat(res.nd_pre_s(:));
% sh_by_unsh_corr = cell2mat(res_corr.nd_pre(:)) ./ cell2mat(res_corr.nd_pre_s(:));
% sh_by_unsh_rnd = cell2mat(res.nd_rnd(:)) ./ cell2mat(res.nd_rnd_s(:));
% sh_by_unsh_rnd_corr = cell2mat(res_corr.nd_rnd(:)) ./ cell2mat(res_corr.nd_rnd_s(:));
% boxplot([sh_by_unsh, sh_by_unsh_corr, sh_by_unsh_rnd(1:38), sh_by_unsh_rnd_corr(1:38)],...
%     {'Unshuffled/Shuffled', 'Unshuffled/Shuffled (corr)', 'Unshuffled/Shuffled Random', 'Unshuffled/Shuffled Random (corr)'});

%%
res = DecodeTensor.noise_properties(T, d, true);
[total_neurons, n_bins, n_trials] = size(T);
%%
figure;%TODO add errorbars/shadedErrorBars from SEM of mean op
s = cell2mat(res.el_pre(:)).^2;
m = mean(s);
e = std(s) ./ sqrt(size(s,1));
shadedErrorBar(1:total_neurons, m, e.*norminv((1+0.95)/2), 'lineprops', 'b');
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');
hold on;
s = cell2mat(res.el_pre_s(:)).^2;
m = mean(s);
e = std(s) ./ sqrt(size(s,1));
shadedErrorBar(1:total_neurons, m, e.*norminv((1+0.95)/2), 'lineprops', 'r');
xlabel 'Noise PC index'
ylabel(sprintf('Signal direction\nin PC basis\n(squared)'))
l_ = refline(0, 1./total_neurons);
l_.Color = 'k';
xlim([1 total_neurons]);
figure_format;
print_svg('signal_pc_loadings');

figure;
n = 1:n_trials/2;
s = cell2mat(res.noise_spectrum(:).').';
m = mean(s); m = m(n);
e = std(s) ./ sqrt(size(s,1)); e = e(n);
shadedErrorBar(n, m, e.*norminv((1+0.95)/2), 'lineprops', 'b');
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');
hold on;
n = 1:n_trials/2;
s = cell2mat(res.noise_spectrum_s(:).').';
m = mean(s); m = m(n);
e = std(s) ./ sqrt(size(s,1)); e = e(n);
shadedErrorBar(n, m, e.*norminv((1+0.95)/2), 'lineprops', 'r');
ylim([1e-2 Inf]);
xlabel 'Noise PC index'
ylabel(sprintf('Correlation matrix\neigenvalue'))
figure_format;
print_svg('noise_spectrum');

figure;
sh_by_unsh = cell2mat(res.nd_pre(:)) ./ cell2mat(res.nd_pre_s(:));
sh_by_unsh_rnd = cell2mat(res.nd_rnd(:)) ./ cell2mat(res.nd_rnd_s(:));
boxplot([sh_by_unsh, sh_by_unsh_rnd(1:38)],...
    {'Signal', 'Random'}, 'outliersize', 1);
l_ = refline(0, 1);
l_.Color = 'k';
ylim([0.3 Inf]);
set(gca, 'YTick', [1 10])
set(gca, 'YScale', 'log');
ylabel(sprintf('Unshuf./Shuf. noise\nvariance ratio'));
figure_format;
print_svg('noise_attenuation_boxplot');