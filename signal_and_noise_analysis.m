recompute = false;

load adjacent_metrics_agg_191202-163306_0.mat

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
res_lookup = @(code) arrayfun(@(x)cellfun(@(y)median(y{code(1),code(2:3)}), x.results_table),res,'UniformOutput',false);
r_ = res_lookup;
medify = @(z) Utils.cf_(@(y)cellfun(@(x)median(x(:)), y), z);
cutoff = 100;
asymp_line = @(n,m) Utils.fitaline(n,m,cutoff);
intercept = @(n,m) Utils.fitaline(n,m,0,true);
show_mice = {'Mouse2022'};%, 'Mouse2024', 'Mouse2028'};
%[~, m_, sp_] = DecodeTensor.special_sess_id_list;
%show_mice = m_;
filt_sess_indices = select_from_mice(show_mice);

%middle
mouse_names = Utils.cf_(@(x)x(17:25), {res.source});
n = {res.n_sizes};

save_fname = 'summarized_adjacent_metrics.mat';
if ~exist(save_fname, 'file')
    c1 = 'fdm';
    c2 = 'scu';
    c3 = 'er';
    progressbar('decoder', 'metric', 'tr/te');
    for i1 = 1:numel(c1)
        c1_i = c1(i1);
        for i2 = 1:numel(c2)
            c2_i = c2(i2);
            for i3 = 1:numel(c3)
                c3_i = c3(i3);
                code = [c1_i c2_i c3_i];
                summarized_adjacent_metrics.(code) = res_lookup(code);
                progressbar(i1/numel(c1), i2/numel(c2), i3/numel(c3));
            end
        end
    end
    save(save_fname, 'summarized_adjacent_metrics');
else
    load(save_fname);
end
s_ = summarized_adjacent_metrics;
%% checking all options:
filt_sess_indices = select_from_mice({'Mouse2022'});

p = Pub(9, 14, 'rows', 3, 'columns', 2);
%dff2_lim = [0 3];

p.panel(1, 'xlab', 'Number of cells', 'ylab', 'Signal ([\Delta{\itF}/{\itF}]^2)', 'title', '(\Delta{\it\mu})^2, on test set');
MultiSessionVisualizer.plot_single_filtered(n, {s_.mse, s_.dse, s_.fse},...
    {'k','m','b'}, filt_sess_indices);
text(150, 1, 'along \Delta{\it\mu} direction', 'Color', 'k', 'HorizontalAlignment', 'center');
text(130, 0.8, '$\hat{w_d}$', 'Color', 'm', 'HorizontalAlignment', 'center', 'Interpreter', 'latex');
text(130, 0.6, '$\hat{w}$', 'Color', 'b', 'HorizontalAlignment', 'center', 'Interpreter', 'latex');
%ylim(dff2_lim);

p.panel(2, 'xlab', 'Number of cells', 'ylab', 'Signal ([\Delta{\itF}/{\itF}]^2)', 'title', '(\Delta{\it\mu})^2, on train set');
MultiSessionVisualizer.plot_single_filtered(n, {s_.msr, s_.dsr, s_.fsr},...
    {'k','m','b'}, filt_sess_indices);
%ylim(dff2_lim);

p.panel(3, 'xlab', 'Number of cells', 'ylab', 'Noise ([\Delta{\itF}/{\itF}]^2)', 'title', '{\it\sigma}^2, on test set');
MultiSessionVisualizer.plot_single_filtered(n, {s_.mce, s_.dce, s_.fce},...
    {'k','m','b'}, filt_sess_indices);
%ylim(dff2_lim);

p.panel(4, 'xlab', 'Number of cells', 'ylab', 'Noise ([\Delta{\itF}/{\itF}]^2)', 'title', '{\it\sigma}^2, on train set');
MultiSessionVisualizer.plot_single_filtered(n, {s_.mcr, s_.dcr, s_.fcr},...
    {'k','m','b'}, filt_sess_indices);
%ylim(dff2_lim);

p.panel(5, 'xlab', 'Number of cells', 'ylab', 'Noise ([\Delta{\itF}/{\itF}]^2)', 'title', '{\it\sigma}_{shuf}^2, on test set');
MultiSessionVisualizer.plot_single_filtered(n, {s_.mue, s_.due, s_.fue},...
    {'k','m','b'}, filt_sess_indices);
%ylim(dff2_lim);

p.panel(6, 'xlab', 'Number of cells', 'ylab', 'Noise ([\Delta{\itF}/{\itF}]^2)', 'title', '{\it\sigma}_{shuf}^2, on train set');
MultiSessionVisualizer.plot_single_filtered(n, {s_.mur, s_.dur, s_.fur},...
    {'k','m','b'}, filt_sess_indices);
%ylim(dff2_lim);

p.format;
p.print('supplements_pdf', 'dm_w_wd_comparison');
%%
n_minus_one = Utils.cf_(@(x)x-1, n);

signal = s_.mse;%medify({res.m2_d});
[signal_slope, signal_slope_conf] = cellfun(asymp_line, n, signal, 'UniformOutput', false);
[signal_intercept, signal_intercept_conf] = cellfun(intercept, n_minus_one, signal, 'UniformOutput', false);

noise = s_.mce;%medify({res.nv_d});
[noise_slope, noise_slope_conf] = cellfun(asymp_line, n, noise, 'UniformOutput', false);
[noise_intercept, noise_intercept_conf] = cellfun(intercept, n_minus_one, noise, 'UniformOutput', false);

noise_shuf = s_.mue;%medify({res.nv_s});
[noise_shuf_slope, noise_shuf_slope_conf] = cellfun(asymp_line, n, noise_shuf, 'UniformOutput', false);
clamp = false;
if clamp
    noise_shuf_slope(cellfun(@(n)n < 0,noise_shuf_slope)) = {eps};
end
[noise_shuf_intercept, noise_shuf_intercept_conf] = cellfun(intercept, n_minus_one, noise_shuf, 'UniformOutput', false);

snr = cellfun(@rdivide, signal, noise, 'UniformOutput', false);
snr_shuf = cellfun(@rdivide, signal, noise_shuf, 'UniformOutput', false);

[asymp_snr, asymp_snr_conf] = cellfun(@uncertain_divide, signal_slope, signal_slope_conf, noise_slope, noise_slope_conf);
[asymp_snr_shuf, asymp_snr_shuf_conf] = cellfun(@uncertain_divide, signal_slope, signal_slope_conf, noise_shuf_slope, noise_shuf_slope_conf);

[intercept_snr, intercept_snr_conf] = cellfun(@uncertain_divide, signal_slope, signal_slope_conf, noise_intercept, noise_intercept_conf);
[intercept_snr_shuf, intercept_snr_shuf_conf] = cellfun(@uncertain_divide, signal_slope, signal_slope_conf, noise_shuf_intercept, noise_shuf_intercept_conf);

nsr = cellfun(@rdivide, noise, signal, 'UniformOutput', false);
nsr_shuf = cellfun(@rdivide, noise_shuf, signal, 'UniformOutput', false);

[asymp_nsr, asymp_nsr_conf] = cellfun(@uncertain_divide, noise_slope, noise_slope_conf, signal_slope, signal_slope_conf);
[asymp_nsr_shuf, asymp_nsr_shuf_conf] = cellfun(@uncertain_divide, noise_shuf_slope, noise_shuf_slope_conf, signal_slope, signal_slope_conf);

% % signal = s_.mse;
% % noise = s_.mce;
% % noise_shuf = s_.mue;
% % 
% % snr = cellfun(@rdivide, signal, noise, 'UniformOutput', false);
% % snr_shuf = cellfun(@rdivide, signal, noise_shuf, 'UniformOutput', false);
%% plotting
filt_sess_indices = select_from_mice({'Mouse2022', 'Mouse2019', 'Mouse2026'});%, 'Mouse2024', 'Mouse2028'});

%p = Pub(12, 7, 'rows', 2, 'columns', 3);
p = Pub(14, 8, 'rows', 2, 'columns', 3);
dff2_lim = [0 1.5];

p.panel(1, 'xlab', 'Number of cells', 'ylab', 'Signal ([\Delta{\itF}/{\itF}]^2)');
MultiSessionVisualizer.plot_single_filtered(n, {signal}, {'k'}, filt_sess_indices);
ylim(dff2_lim);

p.panel(2, 'xlab', 'Number of cells', 'ylab', 'Noise ([\Delta{\itF}/{\itF}]^2)');
MultiSessionVisualizer.plot_single_filtered(n, {noise_shuf, noise}, {'r', 'b'}, filt_sess_indices);
ylim(dff2_lim);
text(100, 1, 'Real', 'Color', 'b');
text(400, 0.4, 'Shuffled', 'Color', 'r');

p.panel(3, 'xlab', 'Number of cells', 'ylab', 'Signal/Noise');
MultiSessionVisualizer.plot_single_filtered(n, {snr, snr_shuf}, {'b', 'r'}, filt_sess_indices);
text(250, 6.2, 'Shuffled', 'Color', 'r');
text(400, 2, 'Real', 'Color', 'b');

y_sh = 0.135;
%p.panel(4, 'ylab', 'Signal slope / Noise slope', 'y_shift', y_sh);

%p.panel(4, 'xlab', 'Mouse index', 'ylab', 'Noise slope / Signal slope', 'y_shift', y_sh);
p.panel(4, 'xlab', 'Mouse index', 'ylab', 'Signal slope / Noise slope', 'y_shift', y_sh);
%[y_,e_] = Utils.bns_groupings(asymp_snr, asymp_snr_shuf, asymp_snr_conf, asymp_snr_shuf_conf, mouse_names, true, {'Real', 'Shuffled'}, true);
%[y_,e_] = Utils.bns_groupings(asymp_nsr, asymp_nsr_shuf, asymp_nsr_conf, asymp_nsr_shuf_conf, mouse_names, false, {'Real', 'Shuffled'}, false, false);
% TODO: CONVERT PANEL 4 TO BOXPLOT OR VIOLIN PLOT OVERLAYED REAL WITH SHUF

hold on;
bf_ = @(x, c) boxplot(x, Utils.cf_(@(x)x(end-1:end),mouse_names), 'Colors', c, 'BoxStyle', 'outline', 'Symbol', '');
%bf_(asymp_nsr, 'b');
bf_(asymp_snr, 'b');
hold on;
l_ = refline(0, 0); l_.Color = 'k'; l_.LineStyle = '-';
%bf_(asymp_nsr_shuf, 'r');
bf_(asymp_snr_shuf, 'r');
hold on;

%ylim([-Inf 6]);
xtickangle(45);
%%%later addition
set(gca, 'YScale', 'log');
%%%


p.panel(5, 'xlab', 'Single cell signal / Single cell noise',...
    'ylab', sprintf('{\\itI}_0 fit value (cm^{-2}%sneuron^{-1})', Utils.dot), 'y_shift', y_sh);
PanelGenerator.plot_regress_averaged(single_dp2(g_), I0_fit_s(g_),...
    1.96.*single_dp2_sem(g_), I0_conf_s(g_), mouse_names(g_), 'r', 'text_coord', [0.028 0.15e-3]);
xlim([0 Inf]);
ylim([0 1e-3]);
p.format;
Utils.fix_exponent(gca , 'y', 0);

p.panel(6, 'xlab', 'Signal slope / Noise slope', 'ylab', '{\itI}_0{\itN} fit value (cm^{-2})', 'y_shift', y_sh);
InfoLimit = N_fit.*I0_fit;
InfoLimit_conf = abs(InfoLimit).*sqrt((N_conf./N_fit).^2 + (I0_conf./I0_fit).^2);
PanelGenerator.plot_regress_averaged(asymp_snr(g_), InfoLimit(g_),...
    asymp_snr_conf(g_), InfoLimit_conf(g_), mouse_names(g_), 'b',...
    'text_coord', [2.4 0.022]);
xlim([0 Inf]);
ylim([0 0.15]);


p.format;
p.print('figure3_pdf', 'SignalNoise');

%% showing for all sessions
% [sess, mouse_names] = DecodeTensor.filt_sess_id_list;
% PanelGenerator.aux_decoding_curves('supplements_pdf/decoding_curves/multi_asnr.pdf',...
%     sess, mouse_names, n, asymp_nsr, asymp_nsr_shuf,...
%     I0_fit, I0_fit_s, N_fit, N_fit_s, 'b', 'r',...
%     [0 0.16], [2 50], [0 Inf]);

figure;
MultiSessionVisualizer.plot_series(n, {snr_shuf, snr}, {'r', 'b'}, mouse_names, [0 Inf]);
            xlabel 'Number of cells'
            ylabel 'Signal / Noise'
            multi_figure_format;
            Utils.printto('supplements_pdf/decoding_curves', 'multi_snr_curves.pdf');
%%
figure;
Utils.bns_groupings(asymp_nsr, asymp_nsr_shuf, asymp_nsr_conf, asymp_nsr_shuf_conf, mouse_names, false);
ylim([-0.5 12]);
ylabel(sprintf('Noise slope / Signal slope'));
xlabel 'Mouse index'
Utils.specific_format('MBNS');
Utils.printto('supplements_pdf/decoding_curves', 'multi_asymp_nsr.pdf');
%% path analysis: splash zone
dm_asnr = get_asnr(res, s_, 'm', false)'; dm_asnr = zscore(dm_asnr(g_));
w_asnr = get_asnr(res, s_, 'f', false)'; w_asnr = zscore(w_asnr(g_));
wd_asnr = get_asnr(res, s_, 'd', false)'; wd_asnr = zscore(wd_asnr(g_));

ss_snr = zscore(single_dp2(g_)');
I0_val = zscore(I0_fit(g_)');
I0N_val = zscore(InfoLimit(g_)');

%% FIRST PATH ANALYSIS - I0 & I0N
%% interaction between I0 and I0N:
fitlm(I0_val, I0N_val)
%% I0, I0N -> ss_snr
fitlm([I0_val, I0N_val], ss_snr)
%% I0, I0N -> dm_asnr
fitlm([I0_val, I0N_val, ss_snr], dm_asnr)
%% SECOND PATH ANALYSIS - ss_snr & dm_asnr
%% interaction between ss_snr and dm_asnr
fitlm(ss_snr, dm_asnr)
%% ss_snr, dm_asnr -> I0
fitlm([ss_snr, dm_asnr], I0_val)
%% ss_snr, dm_asnr, I0_val -> I0N_val
fitlm([ss_snr, dm_asnr, I0_val], I0N_val)
%%
function [asnr, asnr_conf] = get_asnr(res, s_, code, isshuf)
n = {res.n_sizes};
cutoff = 100;
asymp_line = @(n,m) Utils.fitaline(n,m,cutoff);
signal = s_.([code 'se']);%medify({res.m2_d});
[signal_slope, signal_slope_conf] = cellfun(asymp_line, n, signal, 'UniformOutput', false);
%[signal_intercept, signal_intercept_conf] = cellfun(intercept, n_minus_one, signal, 'UniformOutput', false);

if ~isshuf
    noise = s_.([code 'ce']);%medify({res.nv_d});
    [noise_slope, noise_slope_conf] = cellfun(asymp_line, n, noise, 'UniformOutput', false);
    %[noise_intercept, noise_intercept_conf] = cellfun(intercept, n_minus_one, noise, 'UniformOutput', false);
    [asymp_snr, asymp_snr_conf] = cellfun(@uncertain_divide, signal_slope, signal_slope_conf, noise_slope, noise_slope_conf);
    asnr = asymp_snr;
    asnr_conf = asymp_snr_conf;
else
    noise_shuf = s_.([code 'ue']);%medify({res.nv_s});
    [noise_shuf_slope, noise_shuf_slope_conf] = cellfun(asymp_line, n, noise_shuf, 'UniformOutput', false);
    %[noise_shuf_intercept, noise_shuf_intercept_conf] = cellfun(intercept, n_minus_one, noise_shuf, 'UniformOutput', false);
    [asymp_snr_shuf, asymp_snr_shuf_conf] = cellfun(@uncertain_divide, signal_slope, signal_slope_conf, noise_shuf_slope, noise_shuf_slope_conf);
    asnr = asymp_snr_shuf;
    asnr_conf = asymp_snr_shuf_conf;
end
end


%p1 = @process_one_sess;
function [n, signal, noise, noise_shuf, snr, snr_shuf] = process_one_sess(res)
    n = res.n_sizes;
    signal = cellfun(@(t) median(t{'m','se'}), res.results_table);
    noise = cellfun(@(t) median(t{'m','ce'}), res.results_table);
    noise_shuf = cellfun(@(t) median(t{'m','ue'}), res.results_table);
    snr = signal ./ noise;
    snr_shuf = signal ./ noise_shuf;
end


function filt_indices = select_from_mice(mice_to_show)

[~, m_, sp_] = DecodeTensor.special_sess_id_list;
show_filter = ismember(m_, mice_to_show);
filt_indices = sp_(show_filter);

end


function [quotient, quotient_uncertainty] = uncertain_divide(x, xc, y, yc)
quotient = x./y;
quotient_uncertainty = abs(x./y)*sqrt((xc./x).^2 + (yc./y).^2);
end