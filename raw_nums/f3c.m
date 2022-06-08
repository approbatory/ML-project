function f3c
%% Figure 3c
% Raw numbers for:
% - LineEx1X, LineEx2X, LineEx3X
% - LineEx1Y, LineEx2Y, LineEx3Y
% - ShadedEx1Y, ShadedEx2Y, ShadedEx3Y
rng(0);

data = signal_plot;
make_xlsx(data, 'f3c');
end

%% helper functions
function data = signal_plot
s_inds = SessManager.special_sessions_usable_index({'Mouse2022', 'Mouse2024', 'Mouse2028'});
get_usable_result = @(idx) load(onefile(['records_adjacent_metrics_events_transients/*_' num2str(idx, '%.3d') '_*']));

res = cell(1,numel(s_inds));
for s_ix = 1:numel(s_inds)
    res{s_ix} = get_usable_result(s_inds(s_ix));
end

cutoff = 100;
asymp_line = @(n,m) Utils.fitaline(n,m,cutoff);
res_lookup = @(res, code)...
    cellfun(@(x)...
    cellfun(@(y)...
    median(y{code(1),code(2:3)}),...
    x.results_table),...
    res,'UniformOutput',false);


n = cellfun(@(x)x.n_sizes, res, 'UniformOutput', false);
signal = res_lookup(res, 'mse');
noise = res_lookup(res, 'mce');
noise_shuf = res_lookup(res, 'mue');

snr = cellfun(@rdivide, signal, noise,...
    'UniformOutput', false);
snr_shuf = cellfun(@rdivide, signal, noise_shuf,...
    'UniformOutput', false);

[signal_slope, signal_slope_conf] = cellfun(asymp_line, n, signal,...
    'UniformOutput', false);
[noise_slope, noise_slope_conf] = cellfun(asymp_line, n, noise,...
    'UniformOutput', false);
[noise_shuf_slope, noise_shuf_slope_conf] = cellfun(asymp_line, n, noise_shuf,...
    'UniformOutput', false);

[asymp_snr, asymp_snr_conf] = cellfun(@uncertain_divide, signal_slope, signal_slope_conf,...
    noise_slope, noise_slope_conf);
[asymp_snr_shuf, asymp_snr_shuf_conf] = cellfun(@uncertain_divide, signal_slope, signal_slope_conf,...
    noise_shuf_slope, noise_shuf_slope_conf);

%specific for signal
dff2_lim = [0 24];
figure;
data = my_plot_single_filtered(n, {signal}, {'k'}, [true true true]);
ylim(dff2_lim);
xlabel 'Number of cells'
ylabel 'Signal ([event rate]^2)'
figure_format([1 0.75]);
end

function data = my_plot_single_filtered(n_sizes, series_cell, color_cell, filter_selection)
n_sizes = n_sizes(filter_selection);
series_cell = Utils.cf_(@(x)x(filter_selection), series_cell);
for j = 1:numel(series_cell)
    s_ = series_cell{j};
    c_ = color_cell{j};
    for k = 1:numel(s_)
        d = my_neuseries(n_sizes{k}, s_{k}, c_);
        hold on;
        
        data.("LineEx"+k+"X") = d.LineX;
        data.("LineEx"+k+"Y") = d.LineY;
        data.("ShadedEx"+k+"Y") = d.ShadedY;
    end %session
end %quantity shown
xlim([0 500]);
ylim([0 Inf]);
end%func

function data = my_neuseries(n, s, c)
m = mean(s);
e = std(s) ./ sqrt(size(s,1)) .*norminv((1+0.95)/2);
%e = std(s);
h_ = shadedErrorBar(n, m, e, 'lineprops', c);
%errorbar(n, m, e, c);
data.LineX = n;
data.LineY = m;
data.ShadedY = e;
end