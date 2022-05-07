function f3f
%% Figure 3f
% Depends on file: adjacent_metrics_events_transients_agg_210707-152508_0.mat
% Raw numbers for:
% - LineX, LineY, ShadedBottomY, ShadedTopY
rng(0);

data = signal_by_bin;
make_xlsx(data, 'f3f');
end

%% helper functions
function data = signal_by_bin
t_ = tic;
load('adjacent_metrics_events_transients_agg_210707-152508_0.mat', 'res');
toc(t_)
%%
n_sess = numel(res);

[signal, noise, noise_shuf] = deal(cell(1,n_sess));

n_bin_bounds = 38;

progressbar('sess...');
for sess_i = 1:n_sess
    R = res(sess_i).results_table;
    [n_reps, n_n] = size(R);
    [signal{sess_i}, noise{sess_i}, noise_shuf{sess_i}] = deal(zeros(n_reps, n_n, n_bin_bounds));
    for rep_i = 1:n_reps
        for n_i = 1:n_n
            T = R{rep_i, n_i};
            signal{sess_i}(rep_i, n_i, :) = T{'m', 'se'};
            noise{sess_i}(rep_i, n_i, :) = T{'m', 'ce'};
            noise_shuf{sess_i}(rep_i, n_i, :) = T{'m', 'ue'};
        end
    end
    progressbar(sess_i/n_sess);
end
%%
[s_slope, s_slope_conf, n_slope, n_slope_conf, asnr, asnr_conf] = deal(zeros(numel(res), 38));

progressbar('sessions...');
for ix = 1:numel(res)
    [s_slope(ix,:), s_slope_conf(ix,:), n_slope(ix,:), n_slope_conf(ix,:), asnr(ix,:), asnr_conf(ix,:)] =...
        calculate_slopes_bins(ix, signal, noise, noise_shuf, res);
    progressbar(ix/numel(res));
end
%%
kd = bb2kd((1:38).');
d = kd(:,2);

r_plot = @(a,c) shadedErrorBar(1.5:19.5, median(a),...
    [quantile(a,0.75) ; quantile(a,0.25)], 'lineprops', c);

a = s_slope(:,d==0);
r_plot(a, 'g');

data.LineX = 1.5:19.5;
data.LineY = median(a);
data.ShadedBottomY = quantile(a, 0.25);
data.ShadedTopY = quantile(a, 0.75);

set(gca, 'YScale', 'log');

xlabel 'Spatial bin'
ylabel(sprintf('Signal slope\n(events^2/neuron)'));
end

function snr_data = fetch_snr_data(ix, signal, noise, noise_shuf, res)

s = signal{ix};
x.s_m = squeeze(mean(s)); x.s_s = squeeze(std(s)./sqrt(size(s,1)));

n = noise{ix};
x.n_m = squeeze(mean(n)); x.n_s = squeeze(std(n)./sqrt(size(n,1)));

ns = noise_shuf{ix};
x.ns_m = squeeze(mean(ns)); x.ns_s = squeeze(std(ns)./sqrt(size(ns,1)));


snr = s ./ n;
x.snr_m = squeeze(mean(snr)); x.snr_s = squeeze(std(snr)./sqrt(size(snr,1)));

snrs = s ./ ns;
x.snrs_m = squeeze(mean(snrs)); x.snrs_s = squeeze(std(snrs)./sqrt(size(snrs,1)));

x.nneu = res(ix).n_sizes;

snr_data = x;
snr_data.sess_idx = ix;
end

function kd = bb2kd(bb)
kd = [ceil(bb/2), mod(bb,2)];
end


function old_style_plot(snr_data, n_bin_bounds)
x = snr_data;
nneu = x.nneu;
s_m = x.s_m;
s_s = x.s_s;
n_m = x.n_m;
n_s = x.n_s;
ns_m = x.ns_m;
ns_s = x.ns_s;
snr_m = x.snr_m;
snr_s = x.snr_s;
snrs_m = x.snrs_m;
snrs_s = x.snrs_s;

figure('Units',"inches", 'Position',[3.5603    0.8276   16.2500   16.3448]);

bin_colors = jet(n_bin_bounds/2);

for b = 1:n_bin_bounds
    kd = bb2kd(b);
    k = kd(1);
    switch kd(2)
        case 0
            d = 'right';
            is_right = 1;
        case 1
            d = 'left';
            is_right = 0;
    end
    
    
    subplot(6,2,1*2 + is_right - 1);
    shadedErrorBar(nneu.', s_m(:,b), s_s(:,b).*1.96, 'lineprops', {'Color', bin_colors(k,:)});
    xlabel '# neurons'
    ylabel 'Signal'
    ylim([0 max(s_m(:))]);
    %set(gca, 'YScale', 'log');
    
    subplot(6,2,2*2 + is_right - 1);
    shadedErrorBar(nneu.', n_m(:,b), n_s(:,b).*1.96, 'lineprops', {'Color', bin_colors(k,:)});
    xlabel '# neurons'
    ylabel 'Noise'
    ylim([0 max(n_m(:))]);
    %set(gca, 'YScale', 'log');
    
    subplot(6,2,3*2 + is_right - 1);
    shadedErrorBar(nneu.', ns_m(:,b), ns_s(:,b).*1.96, 'lineprops', {'Color', bin_colors(k,:)});
    xlabel '# neurons'
    ylabel 'Noise (shuffled)'
    ylim([0 max(ns_m(:))]);
    %set(gca, 'YScale', 'log');
    
    subplot(6,2,4*2 + is_right - 1);
    shadedErrorBar(nneu.', snr_m(:,b), snr_s(:,b).*1.96, 'lineprops', {'Color', bin_colors(k,:)});
    xlabel '# neurons'
    ylabel 'SNR'
    ylim([0 max(snr_m(:))]);
    %set(gca, 'YScale', 'log');
    
    subplot(6,2,5*2 + is_right - 1);
    shadedErrorBar(nneu.', snrs_m(:,b), snrs_s(:,b).*1.96, 'lineprops', {'Color', bin_colors(k,:)});
    xlabel '# neurons'
    ylabel 'SNR (shuffled)'
    ylim([0 max(snrs_m(:))]);
    %set(gca, 'YScale', 'log');
end
subplot(6,2,1); title 'Signal left'
subplot(6,2,2); title 'Signal right'

subplot(6,2,3); title 'Noise left'
subplot(6,2,4); title 'Noise right'

subplot(6,2,5); title 'Noise (shuf) left'
subplot(6,2,6); title 'Noise (shuf) right'

subplot(6,2,7); title 'SNR left'
subplot(6,2,8); title 'SNR right'

subplot(6,2,9); title 'SNR (shuf) left'
subplot(6,2,10); title 'SNR (shuf) right'
subplot(6,2,[11 12]);
imagesc(1:19);
axis image
colormap jet
set(gca, 'YTickLabels', []);
xticks(1:19);
set(get(gca, 'XAxis'), 'TickLength', [0 0]);
xlabel 'Spatial bins'
sgtitle(sprintf('Session number %d', x.sess_idx));
end

function osp(ix, signal, noise, noise_shuf, res)
old_style_plot(fetch_snr_data(ix, signal, noise, noise_shuf, res), 38);
end

function [r, c] = c_agg_bins(q)
sd_bins = (1:38).';
kd = bb2kd(sd_bins);
k = kd(:,1);
e_bins = ismember(k, [1 2 18 19]);
fs_bins = ismember(k, [3 4 16 17]);
ms_bins = ismember(k, [5 6 14 15]);
ns_bins = ismember(k, [7 8 12 13]);
c_bins = ismember(k, [9 10 11]);

r(:,1) = mean(q(:,e_bins),2);
r(:,2) = mean(q(:,fs_bins),2);
r(:,3) = mean(q(:,ms_bins),2);
r(:,4) = mean(q(:,ns_bins),2);
r(:,5) = mean(q(:,c_bins),2);

c = {'e', 'fs', 'ms', 'ns', 'c'};
end

function [snr_data, categories] = agg_snr_data(snr_data)
vars = {
    's_m'
    's_s'
    'n_m'
    'n_s'
    'ns_m'
    'ns_s'
    'snr_m'
    'snr_s'
    'snrs_m'
    'snrs_s'
    };

for i = 1:numel(vars)
    my_var = vars{i};
    [snr_data.(my_var), categories] = c_agg_bins(snr_data.(my_var));
end

end

function c_legend(cmap, pos, tt)
%axes('Position', [0.2 0.7 0.6 0.6]);
axes('Position', pos);
imagesc([1 1 2 2 3 3 4 4 5 5 5 4 4 3 3 2 2 1 1]);
colormap(gca, cmap)
axis image
axis off
title(tt);
end

function new_plot_sess(isscaled, ix, signal, noise, noise_shuf, res)
snr_data = fetch_snr_data(ix, signal, noise, noise_shuf, res);
[snr_data, ~] = agg_snr_data(snr_data);
x = snr_data;

c_colors = flipud(winter(5));
c_colors_shuf = flipud(spring(5));

nexttile
for i = 1:5
    hold on;
    if isscaled, sc = x.s_m(end,i); else, sc = 1; end
    shadedErrorBar(x.nneu.', x.s_m(:,i)./sc, x.s_s(:,i)./sc, 'lineprops', {'Color', c_colors(i,:)});
end
xlabel 'Num. neurons'
if isscaled
    ylabel(sprintf('Signal\n(scaled)'));
else
    ylabel(sprintf('Signal\n(events^2)'));
end

nexttile
for i = 1:5
    hold on;
    if isscaled, sc = x.n_m(end,i); else, sc = 1; end
    shadedErrorBar(x.nneu.', x.n_m(:,i)./sc, x.n_s(:,i)./sc, 'lineprops', {'Color', c_colors(i,:)});
    shadedErrorBar(x.nneu.', x.ns_m(:,i)./sc, x.ns_s(:,i)./sc, 'lineprops', {'Color', c_colors_shuf(i,:)});
end
xlabel 'Num. neurons'
if isscaled
    ylabel(sprintf('Noise\n(scaled)'));
else
    ylabel(sprintf('Noise\n(events^2)'));
end

nexttile
for i = 1:5
    hold on;
    if isscaled, sc = x.snr_m(end,i); else, sc = 1; end
    shadedErrorBar(x.nneu.', x.snr_m(:,i)./sc, x.snr_s(:,i)./sc, 'lineprops', {'Color', c_colors(i,:)});
    shadedErrorBar(x.nneu.', x.snrs_m(:,i)./sc, x.snrs_s(:,i)./sc, 'lineprops', {'Color', c_colors_shuf(i,:)});
end
xlabel 'Num. neurons'
if isscaled
    ylabel(sprintf('Signal/Noise\n(scaled)'));
else
    ylabel 'Signal/Noise';
end
end

function [asnr, asnr_conf, asnr_shuf, asnr_shuf_conf, asnr_ratio, asnr_ratio_conf] =...
    calculate_asymp_snr(ix, signal, noise, noise_shuf, res)
snr_data = fetch_snr_data(ix, signal, noise, noise_shuf, res);
[snr_data, ~] = agg_snr_data(snr_data);
x = snr_data;

n_cutoff = 100;
for i = 1:5
    [s_slope, s_slope_conf] = Utils.fitaline(x.nneu, x.s_m(:,i).', n_cutoff);
    [n_slope, n_slope_conf] = Utils.fitaline(x.nneu, x.n_m(:,i).', n_cutoff);
    [ns_slope, ns_slope_conf] = Utils.fitaline(x.nneu, x.ns_m(:,i).', n_cutoff);
    
    [asnr(i), asnr_conf(i)] = uncertain_divide(s_slope, s_slope_conf,...
        n_slope, n_slope_conf);
    [asnr_shuf(i), asnr_shuf_conf(i)] = uncertain_divide(s_slope, s_slope_conf,...
        ns_slope, ns_slope_conf);
    [asnr_ratio(i), asnr_ratio_conf(i)] = uncertain_divide(asnr(i), asnr_conf(i),...
        asnr_shuf(i), asnr_shuf_conf(i));
end
end

function [s_slope, s_slope_conf, n_slope, n_slope_conf, asnr, asnr_conf, out] =...
    calculate_slopes_bins(ix, signal, noise, noise_shuf, res)
snr_data = fetch_snr_data(ix, signal, noise, noise_shuf, res);
%[snr_data, ~] = agg_snr_data(snr_data);
x = snr_data;

n_cutoff = 100;
for i = 1:size(x.s_m,2)
    [s_slope(i), s_slope_conf(i)] = Utils.fitaline(x.nneu, x.s_m(:,i).', n_cutoff);
    [n_slope(i), n_slope_conf(i)] = Utils.fitaline(x.nneu, x.n_m(:,i).', n_cutoff);
    [ns_slope(i), ns_slope_conf(i)] = Utils.fitaline(x.nneu, x.ns_m(:,i).', n_cutoff);
    
    [asnr(i), asnr_conf(i)] = uncertain_divide(s_slope(i), s_slope_conf(i),...
        n_slope(i), n_slope_conf(i));
    [asnr_shuf(i), asnr_shuf_conf(i)] = uncertain_divide(s_slope(i), s_slope_conf(i),...
        ns_slope(i), ns_slope_conf(i));
    [asnr_ratio(i), asnr_ratio_conf(i)] = uncertain_divide(asnr(i), asnr_conf(i),...
        asnr_shuf(i), asnr_shuf_conf(i));
end
out.s_slope = s_slope;
out.s_slope_conf = s_slope_conf;

out.n_slope = n_slope;
out.n_slope_conf = n_slope_conf;

out.asnr = asnr;
out.asnr_conf = asnr_conf;

out.asnr_shuf = asnr_shuf;
out.asnr_shuf_conf = asnr_shuf_conf;

out.asnr_ratio = asnr_ratio;
out.asnr_ratio_conf = asnr_ratio_conf;
end