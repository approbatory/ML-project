%creating figure 1
svg_save_dir = 'figure1_svg';
print_svg = @(name) print('-dsvg', fullfile(svg_save_dir, [name '.svg']));

%% panel A: decorrelation result for Mouse2022/21/24
% Caption: Inverse mean squared error for a neural decoder for place as a
% function of number of cells included in the decoder. Aggregated over 20
% random cell subsets of fixed size, shown for one mouse. 
% Shaded regions represent 95% confidence range.
% Trial shuffled data does not reach a performance saturation at 500 neurons, whereas
% unshuffled data does.
name = 'decorrelation_Mouse2022_24_28';
figure;

dbfile = 'decoding.db';

conn = sqlite(dbfile); %remember to close it
%DecodingPlotGenerator.plot_mice(conn, {'Mouse2022'}, 'shuffled', 'NumNeurons', 'IMSE', 'max');
%DecodingPlotGenerator.plot_mice(conn, {'Mouse2022'}, 'unshuffled', 'NumNeurons', 'IMSE', 'max');
mouse_list = {'Mouse2022', 'Mouse2024', 'Mouse2028'};
for m_i = 1:numel(mouse_list)
    mouse = mouse_list{m_i};
    disp(mouse);
    [n,m,e] = DecodingPlotGenerator.get_errors('NumNeurons', conn, mouse, 'unshuffled', 'IMSE', 'max');
    [ns,ms,es] = DecodingPlotGenerator.get_errors('NumNeurons', conn, mouse, 'shuffled', 'IMSE', 'max');
    hold on;
    DecodingPlotGenerator.errors_plotter(ns,ms,es, 'shuffled');
    DecodingPlotGenerator.errors_plotter(n,m,e, 'unshuffled');
end
xlabel 'Number of cells'
ylabel '1/MSE (cm^{-2})'
xlim([0 500]);
text(200+50, 0.3/5.9, 'Unshuffled', 'Color', 'blue');
text(150, 0.8/5.9, 'Shuffled', 'Color', 'red');
figure_format;

print_svg(name);
%print('-dsvg', fullfile(svg_save_dir, 'A.svg'));
%print('-dpng', '-r900', fullfile(svg_save_dir, 'A.png'));
%body
conn.close;
%% confusion mat %%TODO
for i = 1:12
D_T = DecodeTensor(i, 'rawTraces');

[~, ~, ps, ks, ~] = D_T.basic_decode(false, [], []);
[~, ~, ps_s, ks_s, ~] = D_T.basic_decode(true, [], []);

ps = ceil(ps/2);
ks = ceil(ks/2);
ps_s = ceil(ps_s/2);
ks_s = ceil(ks_s/2);

num_trials(i) = 2*D_T.n_one_dir_trials;
C = confusionmat(ks, ps); %./ (2*D_T.n_one_dir_trials);
C_s = confusionmat(ks_s, ps_s);% ./ (2*D_T.n_one_dir_trials);
C_perfect = confusionmat(ks_s, ks);% ./ (2*D_T.n_one_dir_trials);
CDiff(:,:,i) = C_s - C;
CsC = (C_s - C)./C;
CsC(isnan(CsC)) = 0;
disp(i);
end
%%
figure;
m_CDiff = squeeze(sum(CDiff,3))./sum(num_trials);
imagesc(m_CDiff);
colorbar;
%title('Shuffled - Unshuffled');
colormap(bluewhitered);
xlabel 'Predicted bin'
ylabel 'Correct bin'
axis equal;
xlim([0.5 20.5]); ylim([0.5 20.5]);
figure_format('boxsize', [0.75 0.85]); box on;
print_svg('confusion_diff');
% figure; 
% imagesc(C_perfect); colorbar; colormap(bluewhitered);
% title('Perfect CM');
% xlabel 'Predicted bin'
% ylabel 'Correct bin'
% figure;
% imagesc(confusionmat(ks_s, ps_s));
% colorbar;
% title('Shuffled');
%%  panel: full vs diagonal decoder on one mouse:
name = 'fulldiagonal_Mouse2022_24_28';
figure;

dbfile = 'decoding.db';

conn = sqlite(dbfile); %remember to close it
mouse_list = {'Mouse2022', 'Mouse2024', 'Mouse2028'};
for m_i = 1:numel(mouse_list)
    mouse = mouse_list{m_i};
    [n,m,e] = DecodingPlotGenerator.get_errors('NumNeurons', conn, mouse, 'unshuffled', 'IMSE', 'max');
    [ns,ms,es] = DecodingPlotGenerator.get_errors('NumNeurons', conn, mouse, 'diagonal', 'IMSE', 'max');
    hold on;
    DecodingPlotGenerator.errors_plotter(ns,ms,es, 'diagonal');
    DecodingPlotGenerator.errors_plotter(n,m,e, 'unshuffled');
end
xlabel 'Number of cells'
ylabel '1/MSE (cm^{-2})'
xlim([0 500]);
text(300, 0.18/5.9+0.006, 'Full', 'Color', 'blue');
text(150+130, 0.1/5.9, 'Diagonal', 'Color', 'magenta');
figure_format;


print_svg(name);
conn.close;
%% panel B: decorrelation, pooled on all mice
% Caption: As in <'decorrelation_Mouse2022'>, but aggregated over 8 mice.
% 
% Each mouse has a different total number of neurons recorded.
name = 'decorrelation_pooled';
figure;
dbfile = 'decoding.db';
conn = sqlite(dbfile); %remember to close it
mouse_list = {'Mouse2010', 'Mouse2012', 'Mouse2019', 'Mouse2022',...
                'Mouse2023', 'Mouse2024', 'Mouse2026', 'Mouse2028',...
                'Mouse2025', 'Mouse2011', 'Mouse2029', 'Mouse2021'};
pooled_map = containers.Map('KeyType', 'double', 'ValueType', 'any');
pooled_map_s = containers.Map('KeyType', 'double', 'ValueType', 'any');
for i = 1:numel(mouse_list)
    [n,m,~] = DecodingPlotGenerator.get_errors('NumNeurons', conn, mouse_list{i}, 'unshuffled', 'IMSE', 'max');
    [ns,ms,~] = DecodingPlotGenerator.get_errors('NumNeurons', conn, mouse_list{i}, 'shuffled', 'IMSE', 'max');
    for j = 1:numel(n)
        if pooled_map.isKey(n(j))
            pooled_map(n(j)) = [pooled_map(n(j)) m(j)];
        else
            pooled_map(n(j)) = m(j);
        end
    end
    
    for j = 1:numel(ns)
        if pooled_map_s.isKey(ns(j))
            pooled_map_s(ns(j)) = [pooled_map_s(ns(j)) ms(j)];
        else
            pooled_map_s(ns(j)) = ms(j);
        end
    end
end

neuron_nums = [1 (30:30:400)];
imse_means = arrayfun(@(n) mean(pooled_map(n)), neuron_nums);
imse_errbs = arrayfun(@(n) std(pooled_map(n))./sqrt(length(pooled_map(n))), neuron_nums);
imse_counts = arrayfun(@(n) length(pooled_map(n)), neuron_nums);

imse_means_s = arrayfun(@(n) mean(pooled_map_s(n)), neuron_nums);
imse_errbs_s = arrayfun(@(n) std(pooled_map_s(n))./sqrt(length(pooled_map_s(n))), neuron_nums);
imse_counts_s = arrayfun(@(n) length(pooled_map_s(n)), neuron_nums);
hold on;
DecodingPlotGenerator.errors_plotter(neuron_nums,imse_means_s,imse_errbs_s, 'shuffled');
DecodingPlotGenerator.errors_plotter(neuron_nums,imse_means,imse_errbs, 'unshuffled');
%body
xlabel 'Number of cells'
ylabel '1/MSE (cm^{-2})'
xlim([0 400]);
text(10, 0.75/5.9, 'Unshuffled', 'Color', 'blue');
text(10, 0.9/5.9, 'Shuffled', 'Color', 'red');
figure_format;
print_svg(name);
conn.close;

%% panel: decorrelation, pooled on all mice
% and normed to 1 on the imse perf at 
% n_cells = 180 for shuffled
name = 'decorrelation_pooled_normed';
dbfile = 'decoding.db';
conn = sqlite(dbfile); %remember to close it
%mouse_list = {'Mouse2010', 'Mouse2012', 'Mouse2019', 'Mouse2022',...
%                'Mouse2023', 'Mouse2024', 'Mouse2026', 'Mouse2028',...
%                'Mouse2025', 'Mouse2011', 'Mouse2029', 'Mouse2021'};
mouse_list = {'Mouse2019', 'Mouse2021', 'Mouse2022',...
    'Mouse2028', 'Mouse2025', 'Mouse2024'}; %only mice with >80 trials

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
end
lookup_at_value = @(val, val_col, res_col) cell2mat(...
    cellfun(@(n,m) m(n == val), val_col, res_col,...
    'UniformOutput', false).');

m_at_n = @(n) mean(lookup_at_value(n, nP, mP_normed));
e_at_n = @(n) sqrt(var(lookup_at_value(n, nP, mP_normed)) +...
    mean(lookup_at_value(n, nP, eP_normed).^2));

m_at_n_s = @(n) mean(lookup_at_value(n, nPs, mPs_normed));
e_at_n_s = @(n) sqrt(var(lookup_at_value(n, nPs, mPs_normed)) +...
    mean(lookup_at_value(n, nPs, ePs_normed).^2));

m_at_n_d = @(n) mean(lookup_at_value(n, nPd, mPd_normed));
e_at_n_d = @(n) sqrt(var(lookup_at_value(n, nPd, mPd_normed)) +...
    mean(lookup_at_value(n, nPd, ePd_normed).^2));

figure; hold on;
nn = [1 (30:30:500)];
mm = arrayfun(m_at_n, nn);
mm_s = arrayfun(m_at_n_s, nn);
mm_d = arrayfun(m_at_n_d, nn);
ee = arrayfun(e_at_n, nn);
ee_s = arrayfun(e_at_n_s, nn);
ee_d = arrayfun(e_at_n_d, nn);
%errorbar(nn, mm, ee, 'b');
%errorbar(nn, mm_s, ee_s, 'r');
DecodingPlotGenerator.errors_plotter(nn,mm_s,ee_s, 'shuffled');
DecodingPlotGenerator.errors_plotter(nn,mm,ee, 'unshuffled');
xlabel 'Number of cells'
ylabel(sprintf('1/MSE in\nunits of cells'));
xlim([0 500]);
text(10, 350, 'Unshuffled', 'Color', 'blue');
text(10, 420, 'Shuffled', 'Color', 'red');
%text(10, 490, 'y = x', 'Color', 'black');
%axis equal;
l_ = refline(1,0); l_.Color = 'k';
xlim([0 500]); ylim([0 500]);
figure_format;
print_svg(name);

figure; hold on;
DecodingPlotGenerator.errors_plotter(nn, mm, ee, 'unshuffled');
DecodingPlotGenerator.errors_plotter(nn, mm_d, ee_d, 'diagonal');
xlabel 'Number of cells'
ylabel(sprintf('1/MSE in\nunits of cells'));
l_ = refline(1,0); l_.Color = 'k';
xlim([0 500]); ylim([0 200]);
text(10, 350/2.5, 'Full', 'Color', 'blue');
text(10, 420/2.5, 'Diagonal', 'Color', 'magenta');
%text(10, 490/2.5, 'y = x', 'Color', 'black');
figure_format; name = 'fulldiagonal_pooled_normed';
print_svg(name);
conn.close;
%%
figure;
full_gain_m = cellfun(@(x,y)x-y, mP_normed, mPd_normed, 'UniformOutput', false);
full_gain_e = cellfun(@(x,y)sqrt(x.^2+y.^2), eP_normed, ePd_normed, 'UniformOutput', false);
m_at_n_gain = @(n) mean(lookup_at_value(n, nPd, full_gain_m));
e_at_n_gain = @(n) sqrt(var(lookup_at_value(n, nPd, full_gain_m)) +...
    mean(lookup_at_value(n, nPd, full_gain_e).^2));

DecodingPlotGenerator.errors_plotter(nn, arrayfun(m_at_n_gain, nn), arrayfun(e_at_n_gain, nn),...
    'diff');
xlabel 'Number of cells'
ylabel(sprintf('\\Delta1/MSE\nin units of cells'));
text(10, 100, 'Full - Diagonal', 'Color', 'black');
figure_format; name = 'fulldiagonal_normed_diff';
print_svg(name);
%% panel: full vs diagonal decoder on all mice
name = 'fulldiagonal_pooled';
figure;
dbfile = 'decoding.db';
conn = sqlite(dbfile); %remember to close it
mouse_list = {'Mouse2010','Mouse2012',...
    'Mouse2023','Mouse2026',...
    'Mouse2019','Mouse2028',...
    'Mouse2024','Mouse2022'};
pooled_map = containers.Map('KeyType', 'double', 'ValueType', 'any');
pooled_map_s = containers.Map('KeyType', 'double', 'ValueType', 'any');
for i = 1:numel(mouse_list)
    [n,m,~] = DecodingPlotGenerator.get_errors('NumNeurons', conn, mouse_list{i}, 'unshuffled', 'IMSE', 'max');
    [ns,ms,~] = DecodingPlotGenerator.get_errors('NumNeurons', conn, mouse_list{i}, 'diagonal', 'IMSE', 'max');
    for j = 1:numel(n)
        if pooled_map.isKey(n(j))
            pooled_map(n(j)) = [pooled_map(n(j)) m(j)];
        else
            pooled_map(n(j)) = m(j);
        end
    end
    
    for j = 1:numel(ns)
        if pooled_map_s.isKey(ns(j))
            pooled_map_s(ns(j)) = [pooled_map_s(ns(j)) ms(j)];
        else
            pooled_map_s(ns(j)) = ms(j);
        end
    end
end

neuron_nums = [1 (30:30:400)];
imse_means = arrayfun(@(n) mean(pooled_map(n)), neuron_nums);
imse_errbs = arrayfun(@(n) std(pooled_map(n))./sqrt(length(pooled_map(n))), neuron_nums);
imse_counts = arrayfun(@(n) length(pooled_map(n)), neuron_nums);

imse_means_s = arrayfun(@(n) mean(pooled_map_s(n)), neuron_nums);
imse_errbs_s = arrayfun(@(n) std(pooled_map_s(n))./sqrt(length(pooled_map_s(n))), neuron_nums);
imse_counts_s = arrayfun(@(n) length(pooled_map_s(n)), neuron_nums);
hold on;
DecodingPlotGenerator.errors_plotter(neuron_nums,imse_means_s,imse_errbs_s, 'diagonal');
DecodingPlotGenerator.errors_plotter(neuron_nums,imse_means,imse_errbs, 'unshuffled');
%body
xlabel 'Number of cells'
ylabel '1/MSE (cm^{-2})'
xlim([0 400]);
text(10, 0.4/5.9, 'Full', 'Color', 'blue');
text(10, 0.3/5.9, 'Diagonal', 'Color', 'magenta');
figure_format;
print_svg(name);
conn.close;
%% TODO:: more panels for sample decoding (10 neurons/full, shuf/unshuf)
%%and for fit to In/(1+en)
%% panel c: Fitting the number of cells at which information saturates
% Fitting the inverse MSE of place decoders as a function of number of
% cells to I_0 n / (1 + n/N) with parameters I_0 being the performance gain
% per neuron at small number of neurons, and N being the number of neurons
% at which performance is halved compared with purely linear growth.
% Showing fitted N values for 7 mice, as a function of their total
% available cells, color code as in <'decorrelation_Mouse2022'>. Saturation
% occurs much earlier for unshuffled data. Data points not shown (including
% the 8th mouse) if N value diverged, indicating 1/N = 0. Errorbars are 95%
% confidence intervals.
name = 'saturation_fit';
figure;

dbfile = 'decoding.db';

conn = sqlite(dbfile); %remember to close it
%Took out Mouse2010, it has terrible cells&trials, no estimate of N
mouse_list = {'Mouse2012',...
    'Mouse2023','Mouse2026',...
    'Mouse2019','Mouse2028',...
    'Mouse2024','Mouse2022'};
clear n m ns ms
for m_i = 1:numel(mouse_list)
    [n{m_i},m{m_i},~, data_size(m_i)] = DecodingPlotGenerator.get_errors('NumNeurons', conn, ...
        mouse_list{m_i}, 'unshuffled', 'IMSE', 'max');
    [ns{m_i},ms{m_i},~, data_size_s(m_i)] = DecodingPlotGenerator.get_errors('NumNeurons', conn, ...
        mouse_list{m_i}, 'shuffled', 'IMSE', 'max');

    [fitresult, gof] = createFit_infoSaturation(n{m_i}, m{m_i});
    N_fit_value(m_i) = fitresult.N;
    confidence_intervals = confint(fitresult);
    N_lower(m_i) = confidence_intervals(1,2) - N_fit_value(m_i);
    N_upper(m_i) = confidence_intervals(2,2) - N_fit_value(m_i);
    
    [fitresult_s, gof_s] = createFit_infoSaturation(ns{m_i}, ms{m_i});
    N_fit_value_s(m_i) = fitresult_s.N;
    confidence_intervals_s = confint(fitresult_s);
    N_lower_s(m_i) = confidence_intervals_s(1,2) - N_fit_value_s(m_i);
    N_upper_s(m_i) = confidence_intervals_s(2,2) - N_fit_value_s(m_i);
    
end
%
assert(isequal(n,ns));
total_numbers_of_cells = cellfun(@max, n);
errorbar(total_numbers_of_cells, N_fit_value, N_lower, N_upper, '.b', 'CapSize', 3);
hold on;
lower_bound = N_lower_s + N_fit_value_s;
pos_lb = lower_bound > 0;
%errorbar(total_numbers_of_cells, N_fit_value_s, N_lower_s, N_upper_s, 'or');
errorbar(total_numbers_of_cells(pos_lb), N_fit_value_s(pos_lb), N_lower_s(pos_lb), N_upper_s(pos_lb), '.r', 'CapSize', 3);
ylim([10 2000]);
xlim([150 550]);
xlabel 'Total number of cells'
ylabel 'Fitted N'
set(gca, 'YScale', 'log');
set(gca, 'YTick', [10 100 1000]);
set(gca, 'YMinorTick', 'off');
% subplot(1,2,2);
% errorbar(data_size, N_fit_value, N_lower, N_upper, 'o');
% ylim([0 300]);
% xlim([0 200]);
% xlabel 'Number of trials'
% ylabel 'Fitted N'
figure_format;
print_svg(name);
conn.close;

%% Panel d: The effect of correlations on total neuron's activity
% Using event detected neural traces, the distribution of 
[source_path, mouse_name] = DecodeTensor.default_datasets(4);
opt = DecodeTensor.default_opt;
opt.neural_data_type = 'FST_events';%'spikeDeconv';%'FST_events';
opt.restrict_trials = -1;
opt.first_half = false;
[T, d] = DecodeTensor.tensor_loader(source_path, mouse_name, opt);

[X, ks] = DecodeTensor.tensor2dataset(T~=0, d);
X_shuf = shuffle(X, ks);


figure;
name = 'total_activity_dist';
hold on;
histogram(sum(X_shuf,2), 0:15, 'FaceColor', 'r');
histogram(sum(X,2), 0:15, 'FaceColor', 'b');
text(7, 1500, 'Unshuffled', 'Color', 'blue');
text(7, 1100, 'Shuffled', 'Color', 'red');
xlabel 'Number of active cells'
ylabel 'Frequency'
figure_format;
[h,p] = kstest2(sum(X,2), sum(X_shuf,2));
if h
    fprintf('The total active cells distributions are distinct by KS test, p = %e\n', p);
else
    fprintf('The total active cells distributions have not been shown to be distinct, p = %e\n', p);
end
print_svg(name);
%% Panel d_cdf <-- CDF does not look good
name = 'total_activity_CDF';
figure;
hold on;
%[f, x] = ecdf(sum(X,2));
%[fs, xs] = ecdf(sum(X_shuf,2));
%sub = 1:16;
h2 = cdfplot(sum(X_shuf,2)); h2.Color = 'r';
h = cdfplot(sum(X,2)); h.Color = 'b';

xlabel 'Number of active cells'
ylabel 'Empirical CDF'
text(7, 0.5, 'Unshuffled', 'Color', 'blue');
text(7, 0.3, 'Shuffled', 'Color', 'red');
grid off
title '';
xlim([0 15]);
figure_format;
print_svg(name);
%% Panel e_OLD: effect of correlations on total neuron's activity: pooled over 12 mice
name = 'total_activity_change_pooled';
clear h p
corr_effect = cell(12,1);
num_samples = zeros(12,1);
for m_i = 1:12
    [source_path, mouse_name] = DecodeTensor.default_datasets(m_i);
    opt = DecodeTensor.default_opt; opt.first_half = false;
    opt.neural_data_type = 'FST_events';%'spikeDeconv';%'FST_events';
    opt.restrict_trials = -1;
    [T, d] = DecodeTensor.tensor_loader(source_path, mouse_name, opt);
    
    [X, ks] = DecodeTensor.tensor2dataset(T~=0, d);
    X_shuf = shuffle(X, ks);
    
    num_samples(m_i) = size(X,1);
    [cs, eds] = histcounts(sum(X_shuf,2), 0:30);
    [c,  ed ] = histcounts(sum(X,2), 0:30);
    corr_effect{m_i} = c - cs;
    
%     figure;
%     hold on;
%     histogram(sum(X_shuf,2), 0:15, 'FaceColor', 'r');
%     histogram(sum(X,2), 0:15, 'FaceColor', 'b');
%     text(7, 1500, 'Unshuffled', 'Color', 'blue');
%     text(7, 1100, 'Shuffled', 'Color', 'red');
%     xlabel 'Number of active cells'
%     ylabel 'Frequency'
    
    %figure_format;
    [h(m_i),p(m_i)] = kstest2(sum(X,2), sum(X_shuf,2));
    if h(m_i)
        fprintf('The total active cells distributions are distinct by KS test, p = %e\tNcells=%d\tNtrials=%d\n', p(m_i), size(T,1), size(T,3));
    else
        fprintf('The total active cells distributions have not been shown to be distinct, p = %e\tNcells=%d\tNtrials=%d\n', p(m_i), size(T,1), size(T,3));
    end
end
if all(h==1)
    fprintf('They are all distinct');
    disp(p)
end

figure;
corr_effect = cell2mat(corr_effect)./num_samples;
m_c = mean(corr_effect);
e_c = std(corr_effect) ./ sqrt(size(corr_effect,1));
errorbar(0:15, m_c(1:16), e_c(1:16))
xlabel 'Number of active cells'
ylabel '\Delta Freq. / num. samples'
figure_format;
print_svg(name);

%% Panel e : change in std of active cells dists (+sign test)
name = 'change_in_std_pooled';
clear h p
corr_effect = cell(12,1);
num_samples = zeros(12,1);
num_total_cells = zeros(1,12);
stdev_unshuffled = zeros(1,12);
stdev_shuffled = zeros(1,12);
for m_i = 1:12
    [source_path, mouse_name] = DecodeTensor.default_datasets(m_i);
    opt = DecodeTensor.default_opt; opt.first_half = false;
    opt.neural_data_type = 'FST_events';%'spikeDeconv';%'FST_events';
    opt.restrict_trials = -1;
    [T, d] = DecodeTensor.tensor_loader(source_path, mouse_name, opt);
    
    [X, ks] = DecodeTensor.tensor2dataset(T~=0, d);
    %[X, ks] = DecodeTensor.tensor2dataset(T, d);
    X_shuf = shuffle(X, ks);
    
    num_samples(m_i) = size(X,1);
    num_total_cells(m_i) = size(X,2);
    [cs, eds] = histcounts(sum(X_shuf,2), 0:30);
    [c,  ed ] = histcounts(sum(X,2), 0:30);
    stdev_unshuffled(m_i) = std(sum(X,2));
    stdev_shuffled(m_i) = std(sum(X_shuf,2));
    corr_effect{m_i} = c - cs;
    
    [h(m_i),p(m_i)] = kstest2(sum(X,2), sum(X_shuf,2));
    if h(m_i)
        fprintf('The total active cells distributions are distinct by KS test, p = %e\tNcells=%d\tNtrials=%d\n', p(m_i), size(T,1), size(T,3));
    else
        fprintf('The total active cells distributions have not been shown to be distinct, p = %e\tNcells=%d\tNtrials=%d\n', p(m_i), size(T,1), size(T,3));
    end
end
if all(h==1)
    fprintf('They are all distinct');
    disp(p)
end

p_sign = signtest(stdev_unshuffled, stdev_shuffled, 'Tail', 'right');
%%
figure;
%plot([stdev_unshuffled ; stdev_shuffled]./num_total_cells, '-k.');
if false
    plot([stdev_unshuffled ; stdev_shuffled], '-k.');
    xlim([0.5 2.5]);
    %ylim([3 9]*1e-3);
    ylim([0 4]);
    set(gca, 'XTick', [1 2]);
    set(gca, 'XTickLabel', {'Unshuffled', 'Shuffled'});
    %ylabel(sprintf('\\sigma/total cells\nof num. active cells'));
    ylabel(sprintf('\\sigma of num.\nactive cells'));
    text(1.5, 3.5, '***', 'HorizontalAlignment', 'center');
    figure_format;
end
plot(ones(1,12), stdev_unshuffled - stdev_shuffled, 'k.');
set(gca, 'XTick', 1);
set(gca, 'XTickLabel', {});
ylabel(sprintf('\\Delta\\sigma of num.\nactive cells'));
ylim([0 0.5]);
xlabel('Unshuffled - Shuffled');
%sigstar({{'Unshuffled', 'Shuffled'}}, p_sign);
text(1, 0.475, '***', 'HorizontalAlignment', 'center');
figure_format;
if p_sign < 0.05
    fprintf('Sign test significant, p=%e\n', p_sign);
else
    fprintf('Sign test insignificant, p=%e\n', p_sign);
end
print_svg(name);

%% Panel comparing decoding performance for shuf/unshuf, shown on example trajectory
decode_obj = DecodeTensor(4, 'rawTraces', true); %only use first half of session
[me, mse, ps, ks, model] = decode_obj.basic_decode(false, [], []);
[me_s, mse_s, ps_s, ks_s, model_s] = decode_obj.basic_decode(true, [], []);

%%
o = Analyzer('../linear_track/Mouse2022/Mouse-2022-20150326-linear-track/Mouse-2022-20150326_093722-linear-track-TracesAndEvents.mat');
opt = DecodeTensor.default_opt;
[~,~,tr_start,tr_end,~,~,ks] = DecodeTensor.new_sel(o.data.y.raw.full, opt);
tr_mask = false(numel(ks),1);
for tr_i = 1:numel(tr_start)
    tr_mask(tr_start(tr_i):tr_end(tr_i)) = true;
end
second_half_mask = (1:numel(ks)).' > floor(numel(ks)/2);
fast_and_2nd_half = o.data.mask.fast & second_half_mask;
tr_mask = tr_mask(fast_and_2nd_half);

ks = ks(fast_and_2nd_half);
X = o.data.X.full(fast_and_2nd_half, :);

ps = model.predict(X);
ps_s = model_s.predict(shuffle(X, ks));

%perf
mean_err = mean(abs(ceil(ks(tr_mask)/2) - ceil(ps(tr_mask)/2)));
mean_err_s = mean(abs(ceil(ks(tr_mask)/2) - ceil(ps_s(tr_mask)/2)));
%%
f = figure;
xl_ = [0 3];
%ax1 = subplot(2,1,1);
time = ((1:numel(ks)) - 8392)/20;
plot(time, (ceil(ks/2) - 0.5) * opt.bin_width, '-k'); hold on;
plot(time, (ceil(ps/2) - 0.5) * opt.bin_width, '-b'); 
xlim(xl_);
%ylim([1 20]);
%set(gca, 'XTick', []);
%xlabel 'Frame'; 
%title 'Decoding from unshuffled data'
ylabel 'Position (cm)';
xlabel 'Time (s)';

%ax2 = subplot(2,1,2);
hold on;
plot(time, (ceil(ks/2) - 0.5) * opt.bin_width, '-k'); hold on;
plot(time, (ceil(ps_s/2) - 0.5) * opt.bin_width, '-r'); 
%xlim([1440 1700]);
xlim(xl_);
%ylim([1 20]);
%set(gca, 'XTick', []);
%xlabel 'Frame';
%title 'Decoding from shuffled data'
%ylabel 'Position (cm)';
figure_format('boxsize', [0.8 0.7]*1.05); box on;
% %
% f.Units = 'inches';
% f.Position = [f.Position(1:2), [0.8125 0.585].*0.98*1.6.*[1 2]];
% pause(2);
% ax1.Units = 'inches'; ax2.Units = 'inches';
% ax1.FontSize = 6; ax2.FontSize = 6;
% 
% ax1.LineWidth = 0.5; ax2.LineWidth = 0.5;
% ax1.FontName = 'Helvetica LT Std';
% ax2.FontName = 'Helvetica LT Std';
% ax1.TickLength = [0.02 0.02];
% ax2.TickLength = [0.02 0.02];
% 
print_svg('decode_demo');

%% raw trajectory demo
figure; 
time = ((1:numel(o.data.y.raw.full)) - 26000)/20;
plot(time, 118*normalize(o.data.y.raw.full, 'range'), 'k');
[~, ~, tr_s, tr_e, ~, ~, ~] = DecodeTensor.new_sel(o.data.y.raw.full, DecodeTensor.default_opt);
bad_y = 118*normalize(o.data.y.raw.full, 'range');
for tr_i = 1:numel(tr_s)
    bad_y(tr_s(tr_i):tr_e(tr_i)) = nan;
end
hold on;
plot(time, bad_y, 'r');
xlim([0 4000]/20);
ylabel 'Position (cm)'
xlabel 'Time (s)'
%set(gca, 'XTick', []);
%set(gca, 'YTick', []);
figure_format('boxsize', [0.8125*3 0.585].*0.98);
print_svg('raw_traj');