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

%% precompute single cell dp2
save_fname = 'single_cell_dp2.mat';
if ~exist(save_fname, 'file')
    parfor i = 1:107
        [single_dp2(i), single_dp2_sem(i)] = DecodeTensor.cons_filt(i).single_cell_d_primes2;
        fprintf('*');
    end
    fprintf('\ndone\n');
    save(save_fname, 'single_dp2', 'single_dp2_sem');
else
    load(save_fname);
end


%% splash zone
%prelude
medify = @(z) Utils.cf_(@(y)cellfun(@(x)median(x(:)), y), z);
cutoff = 100;
asymp_line = @(n,m) Utils.fitaline(n,m,cutoff);
intercept = @(n,m) Utils.fitaline(n,m,0,true);
%show_mice = {'Mouse2022', 'Mouse2024', 'Mouse2028'};
[~, m_, sp_] = DecodeTensor.special_sess_id_list;
show_mice = m_;
filt_sess_indices = select_from_mice(show_mice);

%middle
mouse_names = Utils.cf_(@(x)x(17:25), {res.source});
n = {res.n_sizes};
n_minus_one = Utils.cf_(@(x)x-1, n);

signal = medify({res.m2_d});
[signal_slope, signal_slope_conf] = cellfun(asymp_line, n, signal, 'UniformOutput', false);
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


%% plotting
p = Pub(12, 7, 'rows', 2, 'columns', 3);
dff2_lim = [0 1.5];

p.panel(1, 'xlab', 'Number of cells', 'ylab', 'Signal ([\DeltaF/F]^2)');
MultiSessionVisualizer.plot_single_filtered(n, {signal}, {'k'}, filt_sess_indices);
ylim(dff2_lim);

p.panel(2, 'xlab', 'Number of cells', 'ylab', 'Noise ([\DeltaF/F]^2)');
MultiSessionVisualizer.plot_single_filtered(n, {noise, noise_shuf}, {'b', 'r'}, filt_sess_indices);
ylim(dff2_lim);
text(100, 1, 'Real', 'Color', 'b');
text(400, 0.4, 'Shuffled', 'Color', 'r');

p.panel(3, 'xlab', 'Number of cells', 'ylab', 'Signal/Noise');
MultiSessionVisualizer.plot_single_filtered(n, {snr, snr_shuf}, {'b', 'r'}, filt_sess_indices);
text(50, 6.2, 'Shuffled', 'Color', 'r');
text(400, 2, 'Real', 'Color', 'b');

y_sh = 0.135;
p.panel(4, 'ylab', 'Signal slope / Noise slope', 'y_shift', y_sh);
Utils.bns_groupings(asymp_snr, asymp_snr_shuf, asymp_snr_conf, asymp_snr_shuf_conf, mouse_names, true, {'Real', 'Shuffled'}, true);

p.panel(5, 'xlab', 'Single cell signal / Single cell noise',...
    'ylab', sprintf('I_0 fit value (cm^{-2}%sneuron^{-1})', Utils.dot), 'y_shift', y_sh);
PanelGenerator.plot_regress_averaged(single_dp2(g_), I0_fit_s(g_),...
    1.96.*single_dp2_sem(g_), I0_conf_s(g_), mouse_names(g_), 'r', 'text_coord', [0.023 0.05e-3]);
xlim([-Inf Inf]);
ylim([-Inf 1e-3]);
Utils.fix_exponent(gca , 'y', 0);

p.panel(6, 'xlab', 'Signal slope / Noise slope', 'ylab', 'I_0N fit value (cm^{-2})', 'y_shift', y_sh);
InfoLimit = N_fit.*I0_fit;
InfoLimit_conf = abs(InfoLimit).*sqrt((N_conf./N_fit).^2 + (I0_conf./I0_fit).^2);
PanelGenerator.plot_regress_averaged(asymp_snr(g_), InfoLimit(g_),...
    asymp_snr_conf(g_), InfoLimit_conf(g_), mouse_names(g_), 'b',...
    'text_coord', [2 0]);
ylim([-Inf 0.15]);

p.format;
p.print('figure3_pdf', 'SignalNoise');
%%

function filt_indices = select_from_mice(mice_to_show)

[~, m_, sp_] = DecodeTensor.special_sess_id_list;
show_filter = ismember(m_, mice_to_show);
filt_indices = sp_(show_filter);

end


function [quotient, quotient_uncertainty] = uncertain_divide(x, xc, y, yc)
quotient = x./y;
quotient_uncertainty = abs(x./y)*sqrt((xc./x).^2 + (yc./y).^2);
end

