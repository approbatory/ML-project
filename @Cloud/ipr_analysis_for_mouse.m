function [signal_ipr, signal_ipr_sem,...
    spectrum_ipr, spectrum_ipr_sem,...
    asymp_ratio_mean, asymp_ratio_sem] =...
    ipr_analysis_for_mouse(mouse_name, use_corr)

if ~exist('use_corr', 'var')
    use_corr = false;
end
if use_corr
    fprintf('Using corr PCA for mouse %s\n', mouse_name);
else
    fprintf('Using cov PCA for mouse %s\n', mouse_name);
end

[~, m_] = DecodeTensor.filt_sess_id_list;

using_session = strcmp(m_, mouse_name);
assert(any(using_session), 'No such mouse: %s', mouse_name);

[spectrum,...
    signal_variances,...
    med_loadings,...
    med_loadings_s,...
    signal_ipr,...
    signal_ipr_s,...
    spectrum_ipr,...
    spectrum_ipr_s,...
    signal_KL,...
    signal_perplexity,...
    area_cos,...
    area_cos2,...
    dmu_IPR,...
    num_neurons,...
    signal_density] = deal(cell(1, sum(using_session)));
    
idx = 0;
for i = 1:numel(m_)
    if using_session(i)
        idx = idx + 1;
        C = Cloud(i, use_corr);
        [spectrum{idx},...
            signal_variances{idx},...
            med_loadings{idx},...
            med_loadings_s{idx},...
            signal_ipr{idx},...
            signal_ipr_s{idx},...
            spectrum_ipr{idx},...
            spectrum_ipr_s{idx},...
            signal_KL{idx},...
            signal_perplexity{idx},...
            area_cos{idx},...
            area_cos2{idx},...
            dmu_IPR{idx},...
            num_neurons{idx},...
            signal_density{idx}] = ...
            C.pc_signal_variance;
    end
end

[spectrum, spectrum_sem] = summ(spectrum);
[signal_variances, signal_variances_sem] = summ(signal_variances);
[med_loadings, med_loadings_sem] = summ(med_loadings);
[med_loadings_s, med_loadings_s_sem] = summ(med_loadings_s);
[signal_ipr, signal_ipr_sem] = summ(signal_ipr);
[signal_ipr_s, signal_ipr_s_sem] = summ(signal_ipr_s);
[spectrum_ipr, spectrum_ipr_sem] = summ(spectrum_ipr);
[spectrum_ipr_s, spectrum_ipr_s_sem] = summ(spectrum_ipr_s);
[signal_KL, signal_KL_sem] = summ(signal_KL);
[signal_perplexity, signal_perplexity_sem] = summ(signal_perplexity);
[signal_density, signal_density_sem] = summ(signal_density);
%asymp snr values, for comparison
load asymp_snr_values.mat

this_mouse_filter = strcmp(mouse_names, mouse_name);
assert(any(this_mouse_filter));
asymp_snr = asymp_snr(this_mouse_filter);
asymp_snr_shuf = asymp_snr_shuf(this_mouse_filter);
asymp_snr_ratio = asymp_snr ./ asymp_snr_shuf;
asymp_ratio_mean = mean(asymp_snr_ratio);
asymp_ratio_sem = sem(asymp_snr_ratio);

if nargout ~= 0
    return;
end
N_ROWS = 3;
figure('Position', [63,687.666666666667,1574.66666666667,410.333333333333]);
if use_corr
    corcovlab = 'corr. PCA';
else
    corcovlab = 'cov. PCA';
end
suptitle(sprintf('All sessions (n=%d) for %s (%s)', sum(using_session), mouse_name, corcovlab));

subplot(1,N_ROWS,1);
hold on;
errorbar(spectrum, spectrum_sem, '-ob');
errorbar(signal_variances, signal_variances_sem, '-og');
%set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
legend 'PC variance' 'Signal variance from PC'
legend box off
xlabel 'PC index'
ylabel 'Variance (\DeltaF/F)^2'

subplot(1,N_ROWS,2);
hold on;
errorbar(med_loadings, med_loadings_sem', '-ob');
errorbar(med_loadings_s, med_loadings_s_sem, '-or');
legend 'Signal-PC loadings' 'Shuf. loadings'
legend box off
xlabel 'PC index'
ylabel 'cos(\Delta\mu,PC_i)'

subplot(1,N_ROWS,3);
hold on;
errorbar(1:39, signal_density, signal_density_sem, '-ob');
bin_labels = [cellfun(@(x)[num2str(x) 'R'],(num2cell(1:19)), 'UniformOutput', false),...
    '', cellfun(@(x)[num2str(x) 'L'],(num2cell(1:19)), 'UniformOutput', false)];
show_ticks = 1:4:39;
set(gca, 'XTick', show_ticks);
set(gca, 'XTickLabel', bin_labels(show_ticks));
legend 'Signal density'
legend Location best
legend box off
xlabel 'Spatial bin'
ylabel 'Signal density'
return;
subplot(1,N_ROWS,3);
hold on;
errorbar(spectrum_ipr, spectrum_ipr_sem, '-ob');
errorbar(spectrum_ipr_s, spectrum_ipr_s_sem, '-or');
errorbar(signal_ipr, signal_ipr_sem, ':xb');
errorbar(signal_ipr_s, signal_ipr_s_sem', ':xr');
legend 'Variance IPR' 'Var. IPR shuf.' 'Signal IPR' 'Signal IPR shuf.'
legend Location best
legend box off
xlabel 'Bin index'
ylabel 'IPR (# of modes)'
set(gca, 'YScale', 'log');

subplot(1,N_ROWS,4);
hold on;
errorbar(signal_perplexity, signal_perplexity_sem, '-ob');
errorbar(signal_ipr, signal_ipr_sem, '-og');
legend 'Signal perp.' 'Signal IPR'
legend Location best
legend box off
xlabel 'Bin index'
ylabel 'Perp. or IPR (# of modes)'

subplot(1,N_ROWS,5);
hold on;
errorbar(signal_KL, signal_KL_sem, '-ob');
legend 'Signal KL'
legend Location best
legend box off
xlabel 'Bin index'
ylabel 'KL Real vs. Shuf'


end


function [mean_, sem_] = summ(x)
num_samples = numel(x);
num_meas = numel(x{1});
a = zeros(num_samples, num_meas);
for i = 1:num_samples
    a(i,:) = x{i};
end
mean_ = mean(a);
sem_ = std(a) ./ sqrt(num_samples);
end