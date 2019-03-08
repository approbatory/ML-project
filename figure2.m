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

%%
figure;
scatter(XS(:,1), XS(:,2), 4, Utils.colorcode(ceil(ks/2)), 'filled', 'MarkerFaceAlpha', 0.02); hold on;
scatter(origin(1), origin(2), 10, 'k');
xlabel PLS1; ylabel PLS2;
axis equal;
xlim(xl_); ylim(yl_);
set(gca, 'XTick', []);
set(gca, 'YTick', []);
text(-0.03, 0.03, 'Unshuffled', 'Color', 'b')
figure_format;
print_png('PLS_adjacent');

figure;
scatter(XS_s(:,1), XS_s(:,2), 4, Utils.colorcode(ceil(ks/2)), 'filled', 'MarkerFaceAlpha', 0.02); hold on;
scatter(origin_s(1), origin_s(2), 10, 'k');
xlabel PLS1; ylabel PLS2;
axis equal;
xlim(xl_); ylim(yl_);
set(gca, 'XTick', []);
set(gca, 'YTick', []);
text(-0.03, 0.03, 'Shuffled', 'Color', 'r')
figure_format;
print_png('PLS_adjacent_shuf');
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
set(gca, 'YTick', [1e-4, 1./total_neurons, 1e-2]);
set(gca, 'YTickLabel', {'10^{-4}', '1/n', '10^{-2}'});
figure_format;
print_svg('signal_pc_loadings');
%%
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
%% regression of nd_pre with perf
dbfile = 'decoding.db';

conn = sqlite(dbfile); %remember to close it
%
mouse_list = {'Mouse2010', 'Mouse2012', 'Mouse2019', 'Mouse2022',...
                'Mouse2023', 'Mouse2024', 'Mouse2026', 'Mouse2028'};

for i = 1:numel(mouse_list)
    [nP{i}, mP{i}, eP{i}] = ...
        DecodingPlotGenerator.get_errors('NumNeurons', conn,...
        mouse_list{i}, 'unshuffled', 'IMSE', 'max');
    [nPs{i}, mPs{i}, ePs{i}] = ...
        DecodingPlotGenerator.get_errors('NumNeurons', conn,...
        mouse_list{i}, 'shuffled', 'IMSE', 'max');
    [nPd{i}, mPd{i}, ePd{i}] = ...
        DecodingPlotGenerator.get_errors('NumNeurons', conn,...
        mouse_list{i}, 'diagonal', 'IMSE', 'max');
    
    %index_of_180 = find(nP{i} == 180,1);
    %norm_by = mP{i}(index_of_180);%REPLACE WITH LINEAR FIT
    %fit_res = createFit_infoSaturation(nPs{i}, mPs{i});
    %norm_by = fit_res.I_0;
    fit_res = fit(double(nPs{i}), mPs{i}, 'p1*x', 'StartPoint', 0.002);
    norm_by(i) = fit_res.p1;
    mPs_normed{i} = mPs{i} ./ norm_by(i);
    ePs_normed{i} = ePs{i} ./ norm_by(i);
    %fit_res = createFit_infoSaturation(double(nP{i}), mP{i});
    %norm_by(i) = fit_res.I_0;
    %fprintf('s I_0: %f\t us I_0: %f\n', norm_by_s(i), norm_by(i));
    mP_normed{i} = mP{i} ./ norm_by(i);
    eP_normed{i} = eP{i} ./ norm_by(i);
    mPd_normed{i} = mPd{i} ./ norm_by(i);
    ePd_normed{i} = ePd{i} ./ norm_by(i);
    
    %calculating noise props `res` variable
    [source_path, mouse_name] = DecodeTensor.default_datasets(i);
    opt = DecodeTensor.default_opt;
    opt.restrict_trials = -1;
    
    [T, d] = DecodeTensor.tensor_loader(source_path, mouse_name, opt);
    res = DecodeTensor.noise_properties(T, d, true); %true for using_corr
    q__ = cell2mat(res.nd_pre(:)) ./ cell2mat(res.nd_pre_s(:));
    sh_by_unsh(i) = mean(q__);
    sh_by_unsh_errbs(i) = std(q__) ./ sqrt(length(q__));
    RMS_noise_corr(i) = res.RMS_noise_corr;
    RMS_noise_corr_shuf(i) = res.RMS_noise_corr_shuf;
    %^ use this to explain that v
    shuffle_effective_neurons_change(i) = mPs_normed{i}(end) - mP_normed{i}(end);
    shuffle_effective_neurons_change_errbs(i) = sqrt(ePs_normed{i}(end).^2 + eP_normed{i}(end).^2);
end
%%
figure; 
scatter(sh_by_unsh, shuffle_effective_neurons_change, 4, 'k','filled');
[fitresult, gof] = fit(sh_by_unsh.', shuffle_effective_neurons_change.', 'poly1');
hold on; 
plot(fitresult, 'k'); legend off
text(6.5, 150, sprintf('adj. R^2=%.2f', gof.adjrsquare));
xlabel(sprintf('Unshuf./Shuf. noise\nvariance ratio'));
ylabel(sprintf('\\Delta1/MSE\nin units of cells'));
figure_format;
print_svg('noise_improvement_regression');
conn.close;
%%
figure; 
xlabel 'RMS noise correlation ratio'
ylabel(sprintf('\\Delta1/MSE\nin units of cells'));
figure_format;
[xData, yData] = prepareCurveData(RMS_noise_corr./RMS_noise_corr_shuf, shuffle_effective_neurons_change);
%plot(RMS_noise_corr./RMS_noise_corr_shuf, shuffle_effective_neurons_change, 'o')
%text(RMS_noise_corr./RMS_noise_corr_shuf+0.005, shuffle_effective_neurons_change+10, mouse_list)
%refline;
scatter(xData, yData, 4, 'k', 'filled');
[fitresult, gof] = fit(xData, yData, 'poly1');
hold on;
plot(fitresult, 'k'); legend off
text(1.7, 400, sprintf('adj. R^2=%.2f', gof.adjrsquare));
xlabel 'RMS noise correlation ratio'
ylabel(sprintf('\\Delta1/MSE\nin units of cells'));
figure_format;
print_svg('RMS_noise_corr_bad_regression');

%% decoding from PLS vs number of dimensions
o = DecodeTensor(4);
n_reps = 20;
dim_list = [1:10, 15, 20, 30, 50, 100, 200, 300, 400, 497];
clear me mse model me_s mse_s model_s me_p mse_p model_p me_sp mse_sp model_sp
for rep = 1:n_reps
    [me(rep), mse(rep), ~, ~, model{rep}] = o.PLS_decode(-1, false, [], []);
    fprintf('no pls: err: %.2f\n', me(rep));
    [me_s(rep), mse_s(rep), ~, ~, model_s{rep}] = o.PLS_decode(-1, true, [], []);
    fprintf('no pls: err: %.2f (shuf)\n', me_s(rep));
    parfor d_i = 1:numel(dim_list)
        d = dim_list(d_i);
        [me_p(rep,d_i), mse_p(rep,d_i), ~, ~, model_p{rep,d_i}] = o.PLS_decode(d, false, [], []);
        fprintf('pls %dd: err: %.2f\n', d, me_p(rep,d_i));
        [me_sp(rep,d_i), mse_sp(rep,d_i), ~, ~, model_sp{rep,d_i}] = o.PLS_decode(d, true, [], []);
        fprintf('pls %dd: err: %.2f (shuf)\n', d, me_sp(rep,d_i));
    end
end
disp(me);
%% saving results
save decoding_pls_dims.mat dim_list me mse model me_s mse_s model_s me_p mse_p model_p me_sp mse_sp model_sp
%% plotting: only 10 reps
figure;
sum_stats = {@mean, @(x)std(x)./sqrt(size(x,1)).*norminv((1+0.95)/2)};
shadedErrorBar(dim_list, 1./mse_p(1:10,:), sum_stats, 'lineprops', 'b');
shadedErrorBar(dim_list, 1./repmat(mse(1:10)', 1, numel(dim_list)), sum_stats, 'lineprops', 'b:');
hold on;
shadedErrorBar(dim_list, 1./mse_sp(1:10,:), sum_stats, 'lineprops', 'r');
shadedErrorBar(dim_list, 1./repmat(mse_s(1:10)', 1, numel(dim_list)), sum_stats, 'lineprops', 'r:');
set(gca, 'XScale', 'log');
xlabel 'PLS dimension'
ylabel '1/MSE (cm^{-2})'
print_svg('decode_from_PLS');

%%
figure;
sum_stats = {@mean, @(x)std(x)./sqrt(size(x,1)).*norminv((1+0.95)/2)};
shadedErrorBar(dim_list, me_p(1:10,:), sum_stats, 'lineprops', 'b');
shadedErrorBar(dim_list, repmat(me(1:10)', 1, numel(dim_list)), sum_stats, 'lineprops', 'b:');
hold on;
shadedErrorBar(dim_list, me_sp(1:10,:), sum_stats, 'lineprops', 'r');
shadedErrorBar(dim_list, repmat(me_s(1:10)', 1, numel(dim_list)), sum_stats, 'lineprops', 'r:');
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
xlabel 'PLS dimension'
ylabel 'Mean error (cm)'
print_svg('decode_from_PLS_mean_err');

%%
o = DecodeTensor(8);
[my_me, my_mse, ~, ~, my_model] = o.PLS_decode(-1, false, [], []);
fprintf('learned unshuf\n');
%[my_me_p, my_mse_p, ~, ~, my_model_p, stats] = o.PLS_decode(497, false, [], []);
[my_me_s, my_mse_s, ~, ~, my_model_s] = o.PLS_decode(-1, true, [], []);
fprintf('learned shuf\n');
%%
%signal_direction = diff(mean(o.data_tensor(:,:, o.tr_dir==1),3),1,2);
[X, ks] = DecodeTensor.tensor2dataset(o.data_tensor, o.tr_dir);
clear b_bs b_s bs_s r_s
for b = 1:o.opt.n_bins-1
    f_bin = 2*b-1; f_bin_next = f_bin+2;
    learner_index = find((my_model.CodingMatrix(f_bin,:)~=0) &...
        (my_model.CodingMatrix(f_bin_next,:)~=0));
    Beta = my_model.BinaryLearners{learner_index}.Beta;
    Beta_s = my_model_s.BinaryLearners{learner_index}.Beta;
    signal_dir = mean(X(ks==f_bin_next,:)) - mean(X(ks==f_bin,:));
    b_bs(b) = angle_v(Beta, Beta_s, true);
    b_s(b) = angle_v(Beta, signal_dir, true);
    bs_s(b) = angle_v(Beta_s, signal_dir, true);
    r_s(b) = angle_v(randn(size(X,2),1), signal_dir, true);
end
figure;
subplot(2,1,1);
plot(b_bs, 'DisplayName', 'Unshuf vs Shuf');
hold on;
plot(b_s, 'DisplayName', 'Unshuf vs \DeltaSignal');
plot(bs_s, 'DisplayName', 'Shuf vs \DeltaSignal');
plot(r_s, 'DisplayName', 'Random vs \DeltaSignal');
legend Location best
ylim([0 90]);
xlabel 'Position bin'
ylabel(sprintf('Angle between\nclassifiers (degrees)'));
for b = 1:o.opt.n_bins-1
    f_bin = 2*b; f_bin_next = f_bin+2;
    learner_index = find((my_model.CodingMatrix(f_bin,:)~=0) &...
        (my_model.CodingMatrix(f_bin_next,:)~=0));
    Beta = my_model.BinaryLearners{learner_index}.Beta;
    Beta_s = my_model_s.BinaryLearners{learner_index}.Beta;
    signal_dir = mean(X(ks==f_bin_next,:)) - mean(X(ks==f_bin,:));
    b_bs(b) = angle_v(Beta, Beta_s, true);
    b_s(b) = angle_v(Beta, signal_dir, true);
    bs_s(b) = angle_v(Beta_s, signal_dir, true);
    r_s(b) = angle_v(randn(size(X,2),1), signal_dir, true);
end
subplot(2,1,2);
plot(b_bs, 'DisplayName', 'b\_bs');
hold on;
plot(b_s, 'DisplayName', 'b\_s');
plot(bs_s, 'DisplayName', 'bs\_s');
plot(r_s, 'DisplayName', 'r\_s');
%legend;
ylim([0 90]);
xlabel 'Position bin'
ylabel(sprintf('Angle between\nclassifiers (degrees)'));