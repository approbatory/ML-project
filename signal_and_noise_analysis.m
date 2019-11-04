recompute = false;

load adjacent_agg_190725-194911_0.mat

fit_savefile = 'decoding_curves_fits.mat';
if recompute || ~exist(fit_savefile, 'file')
    PanelGenerator.decoding_curves('remake', true, 'recompute', recompute);
end
load(fit_savefile);

good_fit_filter = (I0_conf < 0.5*I0_fit) &...
    (I0_conf_s < 0.5*I0_fit_s) &...
    (N_conf < 0.5*N_fit);% & (N_fit < 500) & (cellfun(@max, {res.n_sizes}) >= 300); %200
g_ = good_fit_filter;

%% splash zone
%prelude
medify = @(z) Utils.cf_(@(y)cellfun(@(x)median(x(:)), y), z);
cutoff = 100;
asymp_line = @(n,m) Utils.fitaline(n,m,cutoff);
intercept = @(n,m) Utils.fitaline(n,m,0,true);
%show_mice = {'Mouse2022', 'Mouse2024', 'Mouse2028'};
show_mice = m_;
filt_sess_indices = select_from_mice(show_mice);

%middle
mouse_names = Utils.cf_(@(x)x(17:25), {res.source});
n = {res.n_sizes};
n_minus_one = Utils.cf_(@(x)x-1, n);

signal = medify({res.m2_d});
[signal_slope, signal_slope_conf] = cellfun(@Utils.fitaline, n, signal, 'UniformOutput', false);
[signal_intercept, signal_intercept_conf] = cellfun(intercept, n_minus_one, signal, 'UniformOutput', false);

noise = medify({res.nv_d});
[noise_slope, noise_slope_conf] = cellfun(asymp_line, n, noise, 'UniformOutput', false);
[noise_intercept, noise_intercept_conf] = cellfun(intercept, n_minus_one, noise, 'UniformOutput', false);

noise_shuf = medify({res.nv_s});
[noise_shuf_slope, noise_shuf_slope_conf] = cellfun(asymp_line, n, noise_shuf, 'UniformOutput', false);
[noise_shuf_intercept, noise_shuf_intercept_conf] = cellfun(intercept, n_minus_one, noise_shuf, 'UniformOutput', false);

snr = cellfun(@rdivide, signal, noise, 'UniformOutput', false);
snr_shuf = cellfun(@rdivide, signal, noise_shuf, 'UniformOutput', false);

[asymp_snr, asymp_snr_conf] = cellfun(@uncertain_divide, signal_slope, signal_slope_conf, noise_slope, noise_slope_conf);
[asymp_snr_shuf, asymp_snr_shuf_conf] = cellfun(@uncertain_divide, signal_slope, signal_slope_conf, noise_shuf_slope, noise_shuf_slope_conf);

[intercept_snr, intercept_snr_conf] = cellfun(@uncertain_divide, signal_slope, signal_slope_conf, noise_intercept, noise_intercept_conf);
[intercept_snr_shuf, intercept_snr_shuf_conf] = cellfun(@uncertain_divide, signal_slope, signal_slope_conf, noise_shuf_intercept, noise_shuf_intercept_conf);

dff2_lim = [0 1.5];

%action
figure(1);
subplot(2,3,1);
MultiSessionVisualizer.plot_single_filtered(n, {signal}, {'k'}, filt_sess_indices);
ylim(dff2_lim);
curve_plot_formatting('DFF^2', 'Signal');

subplot(2,3,2);
MultiSessionVisualizer.plot_single_filtered(n, {noise, noise_shuf}, {'b', 'r'}, filt_sess_indices);
curve_plot_formatting('DFF^2', 'Noise');
ylim(dff2_lim);

subplot(2,3,3);
MultiSessionVisualizer.plot_single_filtered(n, {snr, snr_shuf}, {'b', 'r'}, filt_sess_indices);
curve_plot_formatting('SNR^2', 'Signal/Noise');

subplot(2,3,4);
Utils.bns_groupings(asymp_snr, asymp_snr_shuf, asymp_snr_conf, asymp_snr_shuf_conf, mouse_names, true, {'Real', 'Shuffled'}, true);
ylabel('Signal slope / Noise slope');

%change this to actually measure cell by cell the SNR of that cell
figure(2);
Utils.bns_groupings(intercept_snr, intercept_snr_shuf, intercept_snr_conf, intercept_snr_shuf_conf, mouse_names, false, {'Real', 'Shuffled'}, false);
ylabel('Single cell signal / Single cell noise');

figure(1);
subplot(2,3,5);
PanelGenerator.plot_regress_averaged(intercept_snr_shuf(g_), I0_fit_s(g_),...
    intercept_snr_shuf_conf(g_), I0_conf_s(g_), mouse_names(g_), 'r', 'text_coord', [0.1 0.2e-3]);
xlabel('Single cell signal / Single cell noise');
ylabel('I_0');
xlim([-Inf Inf]);
ylim([-Inf Inf]);

subplot(2,3,6);
InfoLimit = N_fit.*I0_fit;
InfoLimit_conf = abs(InfoLimit).*sqrt((N_conf./N_fit).^2 + (I0_conf./I0_fit).^2);
PanelGenerator.plot_regress_averaged(asymp_snr(g_), InfoLimit(g_),...
    asymp_snr_conf(g_), InfoLimit_conf(g_), mouse_names(g_), 'b',...
    'text_coord', [2 -0.03]);
xlabel('Signal slope / Noise slope');
ylabel('I_0N (cm^{-2})');
%%

function filt_indices = select_from_mice(mice_to_show)

[~, m_, sp_] = DecodeTensor.special_sess_id_list;
show_filter = ismember(m_, mice_to_show);
filt_indices = sp_(show_filter);

end

%maybe use sina's code for matlab plots
function curve_plot_formatting(ylab, titl)
xlabel 'Number of cells'
ylabel(ylab);
title(titl);
%figure_format([1 1.4], 'modfig', false);%([1 1.4]);
end

function [quotient, quotient_uncertainty] = uncertain_divide(x, xc, y, yc)
quotient = x./y;
quotient_uncertainty = abs(x./y)*sqrt((xc./x).^2 + (yc./y).^2);
end

