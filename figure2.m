%creating figure 2
svg_save_dir = 'figure2_svg';
print_svg = @(name) print('-dsvg', fullfile(svg_save_dir, [name '.svg']));
print_png = @(name) print('-dpng', '-r1800', fullfile(svg_save_dir, [name '.png']));
%% PLS representations
m_i = 6;
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
xlim([min(XS(:,1)) max(XS(:,1))]);
xlim([min(XS(:,2)) max(XS(:,2))]);
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
scatter(XS(:,1), XS(:,2), 8, Utils.colorcode(ceil(ks/2)), 'filled', 'MarkerFaceAlpha', 0.02); hold on;
scatter(origin(1), origin(2), 10, 'k');
xlabel PLS1; ylabel PLS2;
axis equal; axis tight;
%xlim(xl_); ylim(yl_);
xl_ = xlim;
yl_ = ylim;
set(gca, 'XTick', []);
set(gca, 'YTick', []);
%text(-0.03, 0.03, 'Unshuffled', 'Color', 'b')
figure_format('boxsize', [1 1.2], 'fontsize', 6);
print_png('PLS_adjacent');

figure;
scatter(XS_s(:,1), XS_s(:,2), 8, Utils.colorcode(ceil(ks/2)), 'filled', 'MarkerFaceAlpha', 0.02); hold on;
scatter(origin_s(1), origin_s(2), 10, 'k');
xlabel PLS1; ylabel PLS2;
axis equal; axis tight;
xlim(xl_); ylim(yl_);
set(gca, 'XTick', []);
set(gca, 'YTick', []);
%text(-0.03, 0.03, 'Shuffled', 'Color', 'r')
figure_format('boxsize', [1 1.2], 'fontsize', 6);
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
s = cell2mat(reshape(res.el_pre',[],1)).^2;
%m = mean(s);
%e = std(s) ./ sqrt(size(s,1));
%shadedErrorBar(1:total_neurons, m, e.*norminv((1+0.95)/2), 'lineprops', 'b');
Utils.neuseries(1:total_neurons, s, 'b');
set(gca, 'XScale', 'log'); %set(gca, 'YScale', 'log');
hold on;
s_s = cell2mat(reshape(res.el_pre_s',[],1)).^2;
%m = mean(s);
%e = std(s) ./ sqrt(size(s,1));
%shadedErrorBar(1:total_neurons, m, e.*norminv((1+0.95)/2), 'lineprops', 'r');
Utils.neuseries(1:total_neurons, s_s, 'r');
xlabel 'Noise PC index'
ylabel(sprintf('Signal direction\nin PC basis\n(squared)'))
l_ = refline(0, 1./total_neurons);
l_.Color = 'k';
xlim([1 total_neurons]);
set(gca, 'YTick', [1e-4, 1./total_neurons, 1e-2]);
set(gca, 'YTickLabel', {'10^{-4}', '1/n', '10^{-2}'});
%figure_format;
%print_svg('signal_pc_loadings');
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

%% growth of signal and noise as a function of number of neurons
dm2 = cell(10,1); sm = dm2; sms = dm2;
progressbar('mouse', 'nsize', 'rep');
for m_i = 1:10
    d = DecodeTensor(m_i);
    n_max = size(d.data_tensor,1);
    n_sizes = [1 (30:30:n_max) n_max];
    n_sizes_save{m_i} = n_sizes;
    n_reps = 20;
    %dm2 = zeros(n_reps, numel(n_sizes)); sm = dm2; sms = dm2;
    for n_i = 1:numel(n_sizes)
        for i = 1:n_reps
            [dm2{m_i}(i, n_i), sm{m_i}(i, n_i), sms{m_i}(i, n_i)] = d.signal_and_noise_descriptors(n_sizes(n_i));
            progressbar([],[],i/n_reps);
        end
        fprintf('Done %d of %d\n', n_sizes(n_i), n_sizes(end));
        progressbar([],n_i/numel(n_sizes));
    end
    progressbar(m_i/10);
end
%%
figure;
for m_i = 6
    Utils.neuseries(n_sizes_save{m_i}, dm2{m_i}, 'k'); hold on;
    Utils.neuseries(n_sizes_save{m_i}, sm{m_i}, 'b');
    Utils.neuseries(n_sizes_save{m_i}, sms{m_i}, 'r');
end
xlabel 'Number of cells'
ylabel(sprintf('Population vector\ndistance'));
%legend '(\Delta\mu)^2' '\sigma^2 along \Delta\mu' '\sigma^2 along \Delta\mu (Shuffled)'
text(20, 2, '(\Delta\mu)^2', 'Color', 'k', 'HorizontalAlignment', 'left');
text(20, 2.5, '\sigma^2 along \Delta\mu', 'Color', 'b', 'HorizontalAlignment', 'left');
text(20, 3, '\sigma^2 along \Delta\mu (Shuffled)', 'Color', 'r', 'HorizontalAlignment', 'left');
figure_format;


print_svg('signal_and_noise_growth');
%% only mus
figure;
for m_i = 1:10
    Utils.neuseries(n_sizes_save{m_i}, dm2{m_i}, 'k'); hold on;
end
xlabel 'Number of cells'
ylabel('(\Delta\mu)^2');
%figure_format;

%% normed to maximal deltamu2 value
figure;
for m_i = [3 6 7 4 5 10]
    [dm2_slope(m_i), dm2_conf(m_i)] = Utils.fitaline(n_sizes_save{m_i}, dm2{m_i}, cutoff_n);
    v = dm2_slope(m_i);%mean(dm2{m_i}(:,end))./n_sizes_save{m_i}(end);
    %v = 1;
    Utils.neuseries(n_sizes_save{m_i}, dm2{m_i}./v, 'k'); hold on;
    Utils.neuseries(n_sizes_save{m_i}, sm{m_i}./v, 'b');
    Utils.neuseries(n_sizes_save{m_i}, sms{m_i}./v, 'r');
end
xlabel 'Number of cells'
ylabel(sprintf('Distance\n(In units of cells)'));
%axis equal
text(20, 370, '(\Delta\mu)^2', 'Color', 'k', 'HorizontalAlignment', 'left');
text(20, 470, '\sigma^2 along \Delta\mu', 'Color', 'b', 'HorizontalAlignment', 'left');
text(20, 570, '\sigma^2 along \Delta\mu (Shuffled)', 'Color', 'r', 'HorizontalAlignment', 'left');
figure_format;
print_svg('signal_and_noise_growth');

%% fitted slopes:
cutoff_n = 100;
for m_i = 1:10
    %v = mean(dm2{m_i}(:,end))./n_sizes_save{m_i}(end);
    [dm2_slope(m_i), dm2_conf(m_i)] = Utils.fitaline(n_sizes_save{m_i}, dm2{m_i}, cutoff_n);
    [sm_slope(m_i), sm_conf(m_i)] = Utils.fitaline(n_sizes_save{m_i}, sm{m_i}./dm2_slope(m_i), cutoff_n);
    [sms_slope(m_i), sms_conf(m_i)] = Utils.fitaline(n_sizes_save{m_i}, sms{m_i}./dm2_slope(m_i), cutoff_n);
    [sms_intercept(m_i), sms_intercept_conf(m_i)] = Utils.fitaline(n_sizes_save{m_i}, sms{m_i}./dm2_slope(m_i), cutoff_n, true);
end

figure;
ballnstick('Unshuffled', 'Shuffled', sm_slope, sms_slope, sm_conf, sms_conf);
line(xlim, [1 1], 'Color', 'k', 'LineStyle', ':');
ylim([-Inf Inf]);
ylabel(sprintf('\\sigma^2 along \\Delta\\mu\nrate of change'));
figure_format;
print_svg('rate_of_change_of_sigma2');
%% TODO: same thing as cell above, but for 107 sessions, from signal_and_noise_final.mat
SnN = load('signal_and_noise_final.mat');
n_sess = length(SnN.dm2_full);
cutoff_n = 100;
dm2_slope = zeros(1,n_sess); dm2_conf = dm2_slope; sm_slope = dm2_slope; sm_conf = dm2_slope; sms_slope = dm2_slope; sms_conf = dm2_slope;
for m_i = 1:n_sess
    n_sizes = SnN.n_sizes_full{m_i};
    dm2 = SnN.dm2_full{m_i};
    sm = SnN.sm_full{m_i};
    sms = SnN.sms_full{m_i};
    [dm2_slope(m_i), dm2_conf(m_i)] = Utils.fitaline(n_sizes, dm2, cutoff_n);
    [sm_slope(m_i), sm_conf(m_i)] = Utils.fitaline(n_sizes, sm./dm2_slope(m_i), cutoff_n);
    [sms_slope(m_i), sms_conf(m_i)] = Utils.fitaline(n_sizes, sms./dm2_slope(m_i), cutoff_n);
end
%%
mouse_ids = cellfun(@(x)x(17:25), sheet_paths, 'UniformOutput', false);
mouse_ids = mouse_ids(sheet_paths_filt);
high_conf_filt = (sm_conf < 0.1) & (sms_conf < 0.1); h_ = high_conf_filt;
figure;
ballnstick('Unshuffled', 'Shuffled', sm_slope(h_), sms_slope(h_),...
    sm_conf(h_), sms_conf(h_), 'scatter', true,...
    'coloring', mouse_ids(h_), 'err_scale', 0.05);
line(xlim, [1 1], 'Color', 'k', 'LineStyle', ':');
ylim([-Inf Inf]);
ylabel(sprintf('\\sigma^2 along \\Delta\\mu\nrate of change'));

figure_format;
print_svg('rate_of_change_of_sigma2_scatter');
%% signal and noise plots
figure;
[~, ~, sp_indices] = DecodeTensor.special_sess_id_list;
for m_i = sp_indices
    [dm2_slope, dm2_conf] = Utils.fitaline(n_sizes_full{m_i}, dm2_full{m_i});
    v = dm2_slope;%mean(dm2{m_i}(:,end))./n_sizes_save{m_i}(end);
    %v = 1;
    if mean(sm_full{m_i}(:,end)./v) > 250
        %fprintf('rogue index = %d\n', m_i);
        continue;
    end
    if n_sizes_full{m_i}(end) < 200
        continue;
    end
    fprintf('Including index %d\n', m_i);
    Utils.neuseries(n_sizes_full{m_i}, dm2_full{m_i}./v, 'k'); hold on;
    Utils.neuseries(n_sizes_full{m_i}, sm_full{m_i}./v, 'b');
    Utils.neuseries(n_sizes_full{m_i}, sms_full{m_i}./v, 'r');
end
xlabel 'Number of cells'
ylabel(sprintf('Distance\n(In units of cells)'));
%axis equal
text(20, 370, '(\Delta\mu)^2', 'Color', 'k', 'HorizontalAlignment', 'left');
text(20, 470, '\sigma^2 along \Delta\mu', 'Color', 'b', 'HorizontalAlignment', 'left');
text(20, 570, '\sigma^2 along \Delta\mu (Shuffled)', 'Color', 'r', 'HorizontalAlignment', 'left');
figure_format;
print_svg('signal_and_noise_growth');
%% boxplot of final slope fit
figure;
for s_i = 1:numel(n_sizes_full)
    [dm2_slope(s_i), dm2_conf(s_i)] = Utils.fitaline(n_sizes_full{s_i}, dm2_full{s_i});
    v = dm2_slope(s_i);
    %[changer_fit{s_i}, changer_fit_gof{s_i}] = Utils.fit_slopechanger(n_sizes_full{s_i}, mean(sm_full{s_i})./v);
    %[changer_fit_s{s_i}, changer_fit_gof_s{s_i}] = Utils.fit_slopechanger(n_sizes_full{s_i}, mean(sms_full{s_i})./v);
    [sm_slope(s_i), sm_conf(s_i)] = Utils.fitaline(n_sizes_full{s_i}, sm_full{s_i}./v, 100);
    [sms_slope(s_i), sms_conf(s_i)] = Utils.fitaline(n_sizes_full{s_i}, sms_full{s_i}./v, 100);
end
%%
boxplot([sm_slope; sms_slope]', {'Unshuffled';'Shuffled'});
%%
figure;
ballnstick('Unshuffled', 'Shuffled', sm_slope, sms_slope, sm_conf, sms_conf);
%%bookmark abcdefg
%%
figure;
histogram(sm_slope(h_)); hold on;
histogram(sms_slope(h_));
%% median signal direction loadings 
progressbar('Mouse', 'Ensemble Size', 'Samples');
[sess_list, mouse_list] = DecodeTensor.filt_sess_id_list;
n_ = numel(sess_list);
median_loadings = cell(n_,1); median_loadings_s = median_loadings;
for m_i = 1:n_
    d = DecodeTensor.cons_filt(m_i);
    n_max = size(d.data_tensor,1);
    n_sizes = [unique(ceil(10.^(log10(2):0.1:log10(n_max)))) n_max];
    n_sizes_save{m_i} = n_sizes;
    n_reps = 20;
    for n_i = 1:numel(n_sizes)
        for i = 1:n_reps
            [ml, mls] = d.signal_loadings(n_sizes(n_i));
            if length(ml) < 50
                ml = [ml zeros(1, 50 - length(ml))];
                mls = [mls zeros(1, 50 - length(mls))];
            else
                ml = ml(1:50);
                mls = mls(1:50);
            end
            %[median_loadings(i, n_i, :), median_loadings_s(i, n_i, :)] = d.signal_loadings(n_sizes(n_i));
            median_loadings{m_i}(i, n_i, :) = ml;
            median_loadings_s{m_i}(i, n_i, :) = mls;
            progressbar([],[],i/n_reps);
        end
        fprintf('%d:: Done %d of %d\n', m_i, n_sizes(n_i), n_sizes(end));
        progressbar([], n_sizes(n_i)/n_sizes(end));
    end
    progressbar(m_i/n_);
end
%%
figure;
for m_i = 1:10
    mean_median_loadings = squeeze(mean(abs(median_loadings{m_i})));
    mean_median_loadings_s = squeeze(mean(abs(median_loadings_s{m_i})));
    
    %figure;
    subplot(2,10,m_i);
    im_data = log10(mean_median_loadings(2:end-1,1:30));
    padded_im_data = nan(16 - size(im_data,1), size(im_data,2));
    imagesc([im_data;padded_im_data], [-1.6 log10(0.3)]);
    %colorbar
    %title '|cos(PC, \Delta\mu)|'
    
    ylabel 'Ensemble size'
    set(gca, 'FontSize', 6);
    set(gca, 'FontName', 'Helvetica LT Std');
    set(gca, 'TickLength', [0.02 0.02]);
    rectangle('Position',...
        0.5+[0 size(im_data,1) size(im_data,2) (16 - size(im_data,1))],...
        'FaceColor', 'w', 'EdgeColor', 'k', 'LineStyle', 'none');
    set(gca, 'YTickLabel', 30*cellfun(@str2num, get(gca, 'YTickLabel')));
    %figure_format([0.8 1.4]);
    %box on; colorbar off;
    hc = colorbar;
    
    
    %set(gcf, 'Position', [0 0 1.3 1]);
    %print_svg('signal_loadings_by_size');
    box off
    if m_i > 1
        axis off
    end
    
    %figure;
    subplot(2,10,m_i+10);
    im_data = log10(mean_median_loadings_s(2:end-1,1:30));
    imagesc([im_data;padded_im_data], [-1.6 log10(0.3)]);
    set(gca, 'FontSize', 6);
    set(gca, 'FontName', 'Helvetica LT Std');
    set(gca, 'TickLength', [0.02 0.02]);
    xlabel 'Noise PC'
    rectangle('Position',...
        0.5+[0 size(im_data,1) size(im_data,2) (16 - size(im_data,1))],...
        'FaceColor', 'w', 'EdgeColor', 'k', 'LineStyle', 'none');
    set(gca, 'YTickLabel', 30*cellfun(@str2num, get(gca, 'YTickLabel')));
    %colorbar
    %title '|cos(PC, \Delta\mu)|, Shuffled'
    %set(gca, 'XTickLabel', []);
    %set(gca, 'YTickLabel', []);
    %figure_format([0.8 1.4]);
    box off;
    if m_i > 1
        axis off
    end
end
set(gcf, 'Units', 'inches');
set(gcf, 'Position', [8.6146    4.3854    7.5833    3.4896]);
print_svg('signal_loadings_by_size_shuf');
%%
figure;
cutoff_n = 50;
for m_i = 1:10
    %subplot(2,5,m_i);
    mean_median_loadings = squeeze(mean(median_loadings{m_i}.^2));
    mean_median_loadings_s = squeeze(mean(median_loadings_s{m_i}.^2));
    eig_index = zeros(1, numel(n_sizes_save{m_i})); 
    eig_index_s = eig_index;
    for e_i = 1:numel(n_sizes_save{m_i})
        [~, eig_index(e_i)] = max(mean_median_loadings(e_i, :)); 
        [~, eig_index_s(e_i)] = max(mean_median_loadings_s(e_i,:));
    end
    fprintf('m_i = %d: ', m_i); disp(eig_index); fprintf(' ; '); disp(eig_index_s);
    %Utils.neuseries(n_sizes_save{m_i}, median_loadings{m_i}(:,:,eig_index(end)).^2, 'b');
    loadings_by_size = abs(median_loadings{m_i}(:, sub2ind([numel(n_sizes_save{m_i}) 50], 1:numel(n_sizes_save{m_i}), eig_index)));
    [loadings_slope(m_i), loadings_conf(m_i)] = Utils.fitaline(log10(n_sizes_save{m_i}), log10(loadings_by_size), log10(cutoff_n));
    Utils.neuseries(n_sizes_save{m_i}, loadings_by_size, 'b');
    hold on;
    %Utils.neuseries(n_sizes_save{m_i}, median_loadings_s{m_i}(:,:,eig_index(end)).^2, 'r');
    loadings_by_size_s = abs(median_loadings_s{m_i}(:, sub2ind([numel(n_sizes_save{m_i}) 50], 1:numel(n_sizes_save{m_i}), eig_index_s)));
    [loadings_slope_s(m_i), loadings_conf_s(m_i)] = Utils.fitaline(log10(n_sizes_save{m_i}), log10(loadings_by_size_s), log10(cutoff_n));
    Utils.neuseries(n_sizes_save{m_i}, loadings_by_size_s, 'r');
    plot(1:500, sqrt(0.5)./sqrt(1:500), 'k');
    %title(sprintf('m_i = %d', m_i));
    set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'log');
    xlabel 'Number of cells'
    ylabel('max_\alpha |cos(PC_\alpha, \Delta\mu)|');
end
xlim([1 500]);
ylim([-Inf 1]);
figure_format([1 1.4]);
print_svg('max_noise_loading');

%%
figure;
ballnstick('Unsh.', 'Sh.', loadings_slope, loadings_slope_s, loadings_conf, loadings_conf_s);
line(xlim-0.5, [0 0], 'Color', 'k', 'LineStyle', '-');
line(xlim, [-0.5 -0.5], 'Color', 'k', 'LineStyle', ':');
ylim([-Inf Inf]);
%ylabel('Slope past 50 cells');
figure_format([0.8 1]/2);
print_svg('max_noise_loading_inset_slope');

%% decoding imse saturation value
figure;
%scatter(I0_fit_value.*N_fit_value, 1./sm_slope, 4, 'b');
limit_uncertainty = sqrt((I0_fit_value.*(N_upper)).^2 + (N_fit_value.*(I0_upper)).^2);
inv_sm_slope_uncertainty = sm_conf ./ sm_slope.^2;
hold on;
errorbar(I0_fit_value.*N_fit_value, 1./sm_slope, inv_sm_slope_uncertainty, inv_sm_slope_uncertainty,...
    limit_uncertainty, limit_uncertainty, 'LineStyle', 'none', 'Color', 'b', 'CapSize', 1);
[fitresult, adjr2] = Utils.regress_line(I0_fit_value.*N_fit_value, 1./sm_slope);
h_ = plot(fitresult); legend off
h_.Color = 'b';
text(0.03, 8.5, sprintf('adj. R^2 = %.2f', adjr2));
xlabel 'IMSE limit I_0N';
ylabel(sprintf('Inverse \\sigma^2\nrate of change'));
figure_format;
print_svg('imse_limit_regression');
%% just with N
figure;
scatter(N_fit_value, 1./sm_slope, 40, 'b');
hold on;
errorbar(N_fit_value, 1./sm_slope, inv_sm_slope_uncertainty, inv_sm_slope_uncertainty,...
    N_upper, N_upper, 'LineStyle', 'none', 'Color', 'b');
[fitresult, adjr2] = Utils.regress_line(N_fit_value, 1./sm_slope);
plot(fitresult); legend off
text(0.1, 4, sprintf('adj. R^2 = %.2f', adjr2));
xlabel 'N fit value';
ylabel 'Inverse \sigma^2 rate of change';
%% just with I_0
figure;
scatter(I0_fit_value, 1./sm_slope, 40, 'b');
hold on;
errorbar(I0_fit_value, 1./sm_slope, inv_sm_slope_uncertainty, inv_sm_slope_uncertainty,...
    I0_upper, I0_upper, 'LineStyle', 'none', 'Color', 'b');
[fitresult, adjr2] = Utils.regress_line(I0_fit_value, 1./sm_slope);
plot(fitresult); legend off
text(6e-4, 4, sprintf('adj. R^2 = %.2f', adjr2));
xlabel 'I_0 fit value';
ylabel 'Inverse \sigma^2 rate of change';
%% junk
hold on;
errorbar(N_fit_value(f_), sm_slope(f_),...
    sm_conf(f_), sm_conf(f_), N_lower(f_), N_upper(f_),...
    'LineStyle', 'none', 'Color', 'b');

scatter(N_fit_value_s(f_), sms_slope(f_), 40, 'r');
errorbar(N_fit_value_s(f_), sms_slope(f_),...
    sms_conf(f_), sms_conf(f_), N_lower_s(f_), N_upper_s(f_),...
    'LineStyle', 'none', 'Color', 'r');
set(gca, 'XScale', 'log');
xlabel 'N fit value'
ylabel '\sigma^2 rate of change'
%%
figure;
%scatter(I0_fit_value_s, 1./sms_intercept, 'r');
hold on;
inv_intercept_errb = sms_intercept_conf./sms_intercept.^2;
errorbar(I0_fit_value_s, 1./sms_intercept, inv_intercept_errb, inv_intercept_errb, I0_lower_s, I0_upper_s, 'LineStyle', 'none', 'Color', 'r', 'CapSize', 1);
[fitresult, adjr2] = Utils.regress_line(I0_fit_value_s, 1./sms_intercept);
plot(fitresult); legend off
text(1e-4, 0.055, sprintf('adj. R^2 = %.2f', adjr2));
xlabel 'I_0 fit value'
ylabel 'Asymptotic 1/\sigma^2'
set(gca, 'XTickLabel', arrayfun(@Utils.my_fmt, get(gca, 'XTick') ,'UniformOutput', false));
figure_format;
print_svg('I0_value_regression');
%% sampling issues for estimating N and I0
progressbar('d_neu', 'n_reps');
total_n_neu = 500;
for d_neu = 1:50
    for n_reps = 1:1000
        my_I0 = 5e-4;
        my_N  = 100;
        my_n_samp = [1, d_neu:d_neu:total_n_neu, total_n_neu];
        my_true_IMSE = my_I0 .* my_n_samp ./ (1 + my_n_samp ./ my_N);
        noise_level = 0.001*sqrt(20/n_reps);
        my_noisy_IMSE = my_true_IMSE + noise_level*randn(size(my_true_IMSE));
        [my_fitresult, my_gof] = createFit_infoSaturation(my_n_samp, my_noisy_IMSE);
        c_ = confint(my_fitresult);
        rep_err(d_neu, n_reps) = c_(end) - my_fitresult.N;
        n_calc_reps(d_neu, n_reps) = numel(my_n_samp) * n_reps;
        progressbar(d_neu/50, n_reps/1000);
        %figure;
        %errorbar(my_n_samp, my_noisy_IMSE, noise_level*ones(size(my_true_IMSE)));
    end
end
%%
figure;
semilogy(rep_err); 
refline(0, 10);
refline(0, 5);