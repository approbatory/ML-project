%creating figure 1
svg_save_dir = 'figure1_svg';
print_svg = @(name) print('-dsvg', fullfile(svg_save_dir, [name '.svg']));

%% panel A: decorrelation result for Mouse2022
% Caption: Inverse mean squared error for a neural decoder for place as a
% function of number of cells included in the decoder. Aggregated over 20
% random cell subsets of fixed size, shown for one mouse. 
% Shaded regions represent 95% confidence range.
% Trial shuffled data does not reach a performance saturation at 500 neurons, whereas
% unshuffled data does.
name = 'decorrelation_Mouse2022';
figure;

dbfile = 'decoding.db';

conn = sqlite(dbfile); %remember to close it
%DecodingPlotGenerator.plot_mice(conn, {'Mouse2022'}, 'shuffled', 'NumNeurons', 'IMSE', 'max');
%DecodingPlotGenerator.plot_mice(conn, {'Mouse2022'}, 'unshuffled', 'NumNeurons', 'IMSE', 'max');
mouse = 'Mouse2022';
[n,m,e] = DecodingPlotGenerator.get_errors('NumNeurons', conn, mouse, 'unshuffled', 'IMSE', 'max');
[ns,ms,es] = DecodingPlotGenerator.get_errors('NumNeurons', conn, mouse, 'shuffled', 'IMSE', 'max');
hold on;
DecodingPlotGenerator.errors_plotter(ns,ms,es, 'shuffled');
DecodingPlotGenerator.errors_plotter(n,m,e, 'unshuffled');
xlabel 'Number of cells'
ylabel '1/MSE (cm^{-2})'
xlim([0 500]);
text(200, 0.3, 'Unshuffled', 'Color', 'blue');
text(150, 0.8, 'Shuffled', 'Color', 'red');
figure_format;

print_svg(name);
%print('-dsvg', fullfile(svg_save_dir, 'A.svg'));
%print('-dpng', '-r900', fullfile(svg_save_dir, 'A.png'));
%body
conn.close;

%% panel B: decorrelation, pooled on all mice
% Caption: As in <'decorrelation_Mouse2022'>, but aggregated over 8 mice.
% 
% Each mouse has a different total number of neurons recorded.
name = 'decorrelation_pooled';
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
text(10, 0.75, 'Unshuffled', 'Color', 'blue');
text(10, 0.9, 'Shuffled', 'Color', 'red');
figure_format;
print_svg(name);
conn.close;

%% TODO:: more panels for sample decoding (10 neurons/full, shuf/unshuf)
%%and for fit to I/(1+eN)
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
%% Panel e_OLD: effect of correlations on total neuron's activity: pooled over 8 mice
name = 'total_activity_change_pooled';
clear h p
corr_effect = cell(8,1);
num_samples = zeros(8,1);
for m_i = 1:8
    [source_path, mouse_name] = DecodeTensor.default_datasets(m_i);
    opt = DecodeTensor.default_opt;
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
corr_effect = cell(8,1);
num_samples = zeros(8,1);
num_total_cells = zeros(1,8);
stdev_unshuffled = zeros(1,8);
stdev_shuffled = zeros(1,8);
for m_i = 1:8
    [source_path, mouse_name] = DecodeTensor.default_datasets(m_i);
    opt = DecodeTensor.default_opt;
    opt.neural_data_type = 'FST_events';%'spikeDeconv';%'FST_events';
    opt.restrict_trials = -1;
    [T, d] = DecodeTensor.tensor_loader(source_path, mouse_name, opt);
    
    [X, ks] = DecodeTensor.tensor2dataset(T~=0, d);
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
%%
p_sign = signtest(stdev_unshuffled, stdev_shuffled, 'Tail', 'right');
figure;
plot([stdev_unshuffled ; stdev_shuffled]./num_total_cells, '-k.');
xlim([0.5 2.5]);
ylim([3 9]*1e-3);
set(gca, 'XTick', [1 2]);
set(gca, 'XTickLabel', {'Unshuffled', 'Shuffled'});
ylabel(sprintf('\\sigma/total cells\nof num. active cells'));
text(1.5, 8.5e-3, '**', 'HorizontalAlignment', 'center');
figure_format;
%sigstar({{'Unshuffled', 'Shuffled'}}, p_sign);
%figure_format;
if p_sign < 0.05
    fprintf('Sign test significant, p=%e\n', p_sign);
else
    fprintf('Sign test insignificant, p=%e\n', p_sign);
end
print_svg(name);