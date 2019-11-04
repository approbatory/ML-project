[source_path, mouse_name] = DecodeTensor.default_datasets(4);
opt = DecodeTensor.default_opt;
opt.neural_data_type = 'FST_events';%'spikeDeconv';%'FST_events';
opt.restrict_trials = -1;
[T, d] = DecodeTensor.tensor_loader(source_path, mouse_name, opt);
%%
[X, ks] = DecodeTensor.tensor2dataset(T~=0, d);
X_shuf = shuffle(X, ks);

%%
figure;
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
%%
[cs, eds] = histcounts(sum(X_shuf,2), 0:30);
[c,  ed ] = histcounts(sum(X,2), 0:30);

cum_c = cumsum(c)./sum(c);
cum_cs = cumsum(cs)./sum(cs);

figure; plot(cum_c, 'b'); hold on; plot(cum_cs, 'r');
%% KS test - 2sampled

[h,p] = kstest2(sum(X,2), sum(X_shuf,2), 'Tail', 'larger');
%%
figure; plot(c - cs);
%% redoing on all mice
%f1 = figure; hold on;
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
%%
p_sign = signtest(stdev_unshuffled, stdev_shuffled, 'Tail', 'right');
figure;
plot([stdev_unshuffled ; stdev_shuffled]./num_total_cells, '-k.');
xlim([0.5 2.5]);
set(gca, 'XTick', [1 2]);
set(gca, 'XTickLabel', {'Unshuffled', 'Shuffled'});
ylabel(sprintf('\\sigma/total cells\nof num. active cells'));
figure_format;
if p_sign < 0.05
    fprintf('Sign test significant, p=%e\n', p_sign);
else
    fprintf('Sign test insignificant, p=%e\n', p_sign);
end
% figure;
% corr_effect = cell2mat(corr_effect)./num_samples;
% m_c = mean(corr_effect);
% e_c = std(corr_effect) ./ sqrt(size(corr_effect,1));
% errorbar(0:15, m_c(1:16), e_c(1:16))
% xlabel 'Number of active cells'
% ylabel '\Delta Freq. / num. samples'
% figure_format;