function [mean_table, conf_table] = session_properties_dataset(o)

load asymp_snr_values.mat
load single_cell_dp2.mat
load decoding_curves_fits.mat
load performance.mat
load zero_nc_length.mat
load hd_signal_noise_save.mat

[mat, sem_, varnames, mat_agg, sem_agg] = o.predictor_matrix;


%target: asymp_ratio
[asymp_ratio, asymp_ratio_conf] = ...
    uncertain_divide(asymp_snr, asymp_snr_conf,...
    asymp_snr_shuf, asymp_snr_shuf_conf);

%target: inv_asymp_ratio
[inv_asymp_ratio, inv_asymp_ratio_conf] = ...
    uncertain_divide(asymp_snr_shuf, asymp_snr_shuf_conf,...
    asymp_snr, asymp_snr_conf);

%target: hd_asymp_ratio
[hd_asymp_ratio, hd_asymp_ratio_conf] = ...
    uncertain_divide(hd_asymp_snr, hd_asymp_snr_conf,...
    hd_asymp_snr_shuf, hd_asymp_snr_shuf_conf);

%target: hd_inv_asymp_ratio
[hd_inv_asymp_ratio, hd_inv_asymp_ratio_conf] = ...
    uncertain_divide(hd_asymp_snr_shuf, hd_asymp_snr_shuf_conf,...
    hd_asymp_snr, hd_asymp_snr_conf);

%target: invN50 (real)
invN50 = 1./N_fit; invN50_conf = N_conf ./ N_fit.^2;

%target: N50 (real)
N50 = N_fit; N50_conf = N_conf;

%target: invN50_ratio
[invN50_ratio, invN50_ratio_conf] = uncertain_divide(N_fit_s, N_conf_s, N_fit, N_conf);

%target: N50_ratio
[N50_ratio, N50_ratio_conf] = uncertain_divide(N_fit, N_conf, N_fit_s, N_conf_s);

%target: I0N_ratio
[I0N, I0N_c] = uncertain_multiply(I0_fit, I0_conf, N_fit, N_conf);
[I0N_s, I0N_c_s] = uncertain_multiply(I0_fit_s, I0_conf_s, N_fit_s, N_conf_s);
[I0N_ratio, I0N_ratio_conf] = uncertain_divide(I0N, I0N_c, I0N_s, I0N_c_s);

log_I0N = log10(I0N);
log_I0N_c = I0N_c ./ I0N ./ log(10);

%target: IMSE_150_ratio
chosen_size = 150;
for i = 1:numel(sess)
    ix = find(n_sizes{i} == chosen_size, 1);
    assert(~isempty(ix));
    imse_at_size(i) = mean(imse{i}(:,ix));
    imse_at_size_conf(i) = sem(imse{i}(:,ix));
    
    imse_s_at_size(i) = mean(imse_s{i}(:,ix));
    imse_s_at_size_conf(i) = sem(imse_s{i}(:,ix));
end
[IMSE_150_ratio, IMSE_150_ratio_conf] = ...
    uncertain_divide(imse_at_size, imse_at_size_conf,...
    imse_s_at_size, imse_s_at_size_conf);
IMSE_150_ratio_conf = 1.96 * IMSE_150_ratio_conf;

single_dp2_conf = 1.96 * single_dp2_sem;

%different definitions of performance:
perf_turnaround_only = (number_of_correct_trials + number_of_slow_trials) ./ (number_of_correct_trials + number_of_slow_trials + number_of_turnaround_trials);
perf_slow_only = (number_of_correct_trials + number_of_turnaround_trials) ./ (number_of_correct_trials + number_of_slow_trials + number_of_turnaround_trials);
%make tables from the variables:


mean_table = array2table([mat, asymp_ratio', inv_asymp_ratio', invN50', N50',...
    invN50_ratio', N50_ratio', I0N_ratio', IMSE_150_ratio', I0N', asymp_snr', N_fit_s', I0_fit', I0_fit_s', I0N_s', asymp_snr_shuf', single_dp2', performance', perf_turnaround_only', perf_slow_only', zero_nc_length', log_I0N', hd_asymp_ratio, hd_inv_asymp_ratio],...
    'VariableNames', [varnames, {'asymp_ratio', 'inv_asymp_ratio', 'invN50', 'N50',...
    'invN50_ratio', 'N50_ratio', 'I0N_ratio', 'IMSE_150_ratio', 'I0N', 'asymp_snr', 'N50_s', 'I0', 'I0_s', 'I0N_s', 'asymp_snr_shuf', 'single_dp2', 'performance', 'perf_turnaround_only', 'perf_slow_only', 'zero_nc_length', 'log_I0N', 'hd_asymp_ratio', 'hd_inv_asymp_ratio'}]);
mean_table = [table(mouse_names', 'VariableNames', {'Mouse'}), mean_table];

conf_table = array2table([1.96*sem_, asymp_ratio_conf', inv_asymp_ratio_conf', invN50_conf', N50_conf',...
    invN50_ratio_conf', N50_ratio_conf', I0N_ratio_conf', IMSE_150_ratio_conf', I0N_c', asymp_snr_conf', N_conf_s', I0_conf', I0_conf_s', I0N_c_s', asymp_snr_shuf_conf', single_dp2_conf', performance'.*0, perf_turnaround_only'.*0, perf_slow_only'.*0, zero_nc_length'.*0, log_I0N_c.', hd_asymp_ratio_conf, hd_inv_asymp_ratio_conf],...
    'VariableNames', [varnames, {'asymp_ratio', 'inv_asymp_ratio', 'invN50', 'N50',...
    'invN50_ratio', 'N50_ratio', 'I0N_ratio', 'IMSE_150_ratio', 'I0N', 'asymp_snr', 'N50_s', 'I0', 'I0_s', 'I0N_s', 'asymp_snr_shuf', 'single_dp2', 'performance', 'perf_turnaround_only', 'perf_slow_only', 'zero_nc_length', 'log_I0N', 'hd_asymp_ratio', 'hd_inv_asymp_ratio'}]);
conf_table = [table(mouse_names', 'VariableNames', {'Mouse'}), conf_table];


o.sess_prop = mean_table;
o.sess_prop_conf = conf_table;
end



function [quotient, quotient_uncertainty] = uncertain_divide(x, xc, y, yc)
quotient = x./y;
quotient_uncertainty = abs(x./y).*sqrt((xc./x).^2 + (yc./y).^2);
end

function [product, product_uncertainty] = uncertain_multiply(x, xc, y, yc)
product = x.*y;
product_uncertainty = abs(product) .* sqrt((xc./x).^2 + (yc./y).^2);
end