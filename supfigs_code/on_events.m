ntype = 'HD'; % 'spikeDeconvTrace', 'IED'

%cd ~/ML-project/
d = DecodeTensor(6, ntype);
%%
p = Pub(14, 16, 'rows', 4, 'columns', 4, 'vspacing', 0.08,...
    'vmargint', 0.05, 'vmarginb', 0.05, 'hspacing', 0.11);
%p.preview;
%%
alpha = 0.02;

[X, ks] = d.get_dataset;
X_shuf = shuffle(X, ks);
X_z = zscore(X);
X_shuf_z = zscore(X_shuf);
[XS, stats, origin] = Utils.pls_short(X_z, [ceil(ks/2), mod(ks,2)]);
XS_s = (X_shuf_z - mean(X_shuf_z)) * stats.W;
origin_s = -mean(X_shuf_z) * stats.W;

p.panel(1, 'xlab', 'PLS1', 'ylab', 'PLS2', 'y_shift', 0.04);
scatter(XS(:,1), XS(:,2), 8, Utils.colorcode(ceil(ks/2)), 'filled', 'MarkerFaceAlpha', alpha); hold on;
scatter(origin(1), origin(2), 10, 'k');
%xlabel PLS1; ylabel PLS2;
axis equal; axis tight;
xl_ = xlim;
yl_ = ylim;
set(gca, 'XTick', []);
set(gca, 'YTick', []);
%figure_format('boxsize', [1 1.2], 'fontsize', 6);
%print('-dpng', '-r1800', 'figure2_pdf/demo/PLS_adjacent.png');

p.panel(2, 'xlab', 'PLS1', 'ylab', 'PLS2', 'y_shift', 0.04);
scatter(XS_s(:,1), XS_s(:,2), 8, Utils.colorcode(ceil(ks/2)), 'filled', 'MarkerFaceAlpha', alpha); hold on;
scatter(origin_s(1), origin_s(2), 10, 'k');
%xlabel PLS1; ylabel PLS2;
axis equal; axis tight;
xlim(xl_); ylim(yl_);
set(gca, 'XTick', []);
set(gca, 'YTick', []);
%figure_format('boxsize', [1 1.2], 'fontsize', 6);
%print('-dpng', '-r1800', 'figure2_pdf/demo/PLS_adjacent_shuf.png');

% figure;
% my_colors = Utils.colorcode(1:20);
% for b_i = 1:20
%     rectangle('Position', [b_i 0 1 0.5], 'FaceColor', my_colors(b_i,:), 'EdgeColor', 'none');
% end
% axis equal; axis tight; axis off;
p.format;

%%

save_file = sprintf('decoding_curves_fits_%s.mat', ntype);

if ~exist(save_file, 'file')
    dbfile = sprintf('decoding_all_sess_%s.db', ntype);
    conn = sqlite(dbfile);
    samp_size = 80;
    %[sess, mouse_names] = DecodeTensor.filt_sess_id_list;
    [sess, mouse_names] = SessManager.usable_sess_id_list;
    sess = cellfun(@(a,b) {a,b}, mouse_names, sess, 'UniformOutput', false);
    [n_sizes, imse, mask] = PanelGenerator.db_imse_reader_safe(conn, 'unshuffled', sess, samp_size);
    [n_sizes_s, imse_s, mask_s] = PanelGenerator.db_imse_reader_safe(conn, 'shuffled', sess, samp_size);
    [n_sizes_d, imse_d, mask_d] = PanelGenerator.db_imse_reader_safe(conn, 'diagonal', sess, samp_size);
    assert(isequal(n_sizes, n_sizes_s), 'mismatch between unshuffled and shuffled sampling');
    assert(isequal(n_sizes, n_sizes_d), 'mismatch between unshuffled and diagonal sampling');
    
    [series_fits{1}, series_gof{1}] = Utils.cf_p2(2,@(n,m)createFit_infoSaturation(n(:),mean(m)'), n_sizes, imse);
    progressbar(1/3/2);
    [series_fits{2}, series_gof{2}] = Utils.cf_p2(2,@(n,m)createFit_infoSaturation(n(:),mean(m)'), n_sizes, imse_s);
    progressbar(2/3/2);
    [series_fits{3}, series_gof{3}] = Utils.cf_p2(2,@(n,m)createFit_infoSaturation(n(:),mean(m)'), n_sizes, imse_d);
    progressbar(3/3/2);
    [series_exp_fits{1}, series_exp_gof{1}] = Utils.cf_p2(2,@(n,m)createFit_exp(n(:),mean(m)'), n_sizes, imse);
    progressbar(4/3/2);
    [series_exp_fits{2}, series_exp_gof{2}] = Utils.cf_p2(2,@(n,m)createFit_exp(n(:),mean(m)'), n_sizes, imse_s);
    progressbar(5/3/2);
    [series_exp_fits{3}, series_exp_gof{3}] = Utils.cf_p2(2,@(n,m)createFit_exp(n(:),mean(m)'), n_sizes, imse_d);
    progressbar(6/3/2);
    
    [I0_fit, I0_conf] = Utils.fit_get(series_fits{1}, 'I_0');
    [I0_fit_s, I0_conf_s] = Utils.fit_get(series_fits{2}, 'I_0');
    [I0_fit_d, I0_conf_d] = Utils.fit_get(series_fits{3}, 'I_0');
    
    [N_fit, N_conf] = Utils.fit_get(series_fits{1}, 'N');
    [N_fit_s, N_conf_s] = Utils.fit_get(series_fits{2}, 'N');
    [N_fit_d, N_conf_d] = Utils.fit_get(series_fits{3}, 'N');
    
    save(save_file, 'sess', 'mouse_names', 'n_sizes',...
        'imse', 'imse_s', 'imse_d', 'series_fits', 'series_gof', 'I0_fit', 'I0_conf',...
        'I0_fit_s', 'I0_conf_s', 'N_fit', 'N_conf',...
        'N_fit_s', 'N_conf_s', 'series_exp_fits', 'series_exp_gof',...
        'I0_fit_d', 'I0_conf_d', 'N_fit_d', 'N_conf_d', 'mask');
else
    load(save_file);
end

r2 = cellfun(@(x)x.rsquare, series_gof{1});
r2_s = cellfun(@(x)x.rsquare, series_gof{2});
r2_d = cellfun(@(x)x.rsquare, series_gof{3});

r2_exp = cellfun(@(x)x.rsquare, series_exp_gof{1});
r2_s_exp = cellfun(@(x)x.rsquare, series_exp_gof{2});
r2_d_exp = cellfun(@(x)x.rsquare, series_exp_gof{3});

fprintf('Unshuf R^2: %f - %f, median: %f\n', min(r2), max(r2), median(r2));
fprintf('Shuf R^2: %f - %f, median: %f\n', min(r2_s), max(r2_s), median(r2_s));
fprintf('Diag R^2: %f - %f, median: %f\n', min(r2_d), max(r2_d), median(r2_d));
fprintf('Unshuf (exp) R^2: %f - %f, median: %f\n', min(r2_exp), max(r2_exp), median(r2_exp));
fprintf('Shuf (exp) R^2: %f - %f, median: %f\n', min(r2_s_exp), max(r2_s_exp), median(r2_s_exp));
fprintf('Diag (exp) R^2: %f - %f, median: %f\n', min(r2_d_exp), max(r2_d_exp), median(r2_d_exp));


%PanelGenerator.aux_decoding_curves(fname, sess, mouse_names, n_sizes, imse, imse_s,...
%                    I0_fit, I0_fit_s, N_fit, N_fit_s, 'b', 'r',...
%                    [0 0.16], [2 50], [0 Inf]);
%%
selected_indices = SessManager.special_sessions_usable_index({'Mouse2022', 'Mouse2024', 'Mouse2028'});
sm_ = SessManager;
s_inds = false(1,sm_.num_usable); s_inds(selected_indices) = true;
s_inds = s_inds(mask);
s_inds = find(s_inds);

p.panel([3 4], 'xlab', 'Number of cells', 'ylab', '1/MSE (cm^{-2})', 'y_shift', 0.04, 'letter', 'c');
%PanelGenerator.plot_decoding_curve(sess, s_inds, n_sizes, imse_s, I0_fit_s, N_fit_s, 'r');
%hold on;
%PanelGenerator.plot_decoding_curve(sess, s_inds, n_sizes, imse, I0_fit, N_fit, 'b');
MultiSessionVisualizer.plot_single_filtered(n_sizes, {imse_s, imse}, {'r', 'b'}, s_inds);
ylim([0 Inf]);
legend off

%%
%res_collection = arrayfun(@(x)med_loadings_compute(x,'IED'), selected_indices, 'UniformOutput', false);
%save medload_selected_res_collection.mat res_collection selected_indices
%%
% load medload_selected_res_collection.mat
% res = cell2mat(res_collection);
% n_sizes = {res.n_sizes};
% n_sizes = Utils.cf_(@(x)x(2:end),n_sizes);
% series = {{res.median_loadings_s}, {res.median_loadings}};
% 
% series = Utils.cf_(@(m)Utils.cf_(@(x)max(x(:,2:end,:),[],3),m), series);
% mouse_name = {res.mouse_name};
% 
% show_mice = {'Mouse2022', 'Mouse2024', 'Mouse2028'};
% 
% [~,m_,sp_] = DecodeTensor.special_sess_id_list;
% show_filter = ismember(m_, show_mice);
% sp_ = sp_(show_filter);
% m_ = m_(show_filter);
% %figure;
% colorscale = 'log';
% for i = 1:numel(sp_)
%     p.panel(4+i, 'xlab', 'Fluctuation mode, {\iti}', 'ylab', 'Number of cells',...
%         'title', sprintf('Mouse %s', m_{i}(end-1:end)), 'letter', char('c'+i), 'y_shift', 0.04);
%     %subplot(1,numel(sp_)+1, i);
%     mean_median_loadings = squeeze(mean(abs(res(i).median_loadings)));
%     min_d = 30;
%     ns = res(i).n_sizes;
%     im_data = (mean_median_loadings(ns >= min_d,1:min_d));
%     %padded_im_data = nan(16 - size(im_data,1), size(im_data,2));
%     %imagesc([im_data;padded_im_data], [-1.6 log10(0.3)]);
%     surf(1:min_d, ns(ns>=min_d), im_data, 'EdgeColor', 'none');
%     view(2);
%     
%     %set(gca, 'XScale', 'log');
%     %set(gca, 'YScale', 'log');
%     set(gca, 'ColorScale', colorscale);
%     caxis([0.03 0.25]);
%     xlim([1 min_d]);
%     ylim([min_d+10, 500]);
%     %xlabel 'Fluctuation mode, i'
%     %ylabel 'Number of cells'
%     title(sprintf('Mouse %s', m_{i}(end-1:end)), 'FontName', 'Helvetica', 'FontSize', 6, 'FontWeight', 'normal', 'Color', 'b');
%     
%     set(gca, 'FontSize', 6);
%     set(gca, 'FontName', 'Helvetica');
%     set(gca, 'TickLength', [0.02 0.02]);
%     %colorbar;
%     %set(gca, 'YTick', [1 2 4 8 16 32 64].*min_d);
%     %rectangle('Position',...
%     %    0.5+[0 size(im_data,1) size(im_data,2) (48 - size(im_data,1))],...
%     %    'FaceColor', 'w', 'EdgeColor', 'k', 'LineStyle', 'none');
%     %set(gca, 'YTickLabel', 10*cellfun(@str2num, get(gca, 'YTickLabel')));
%     box off;
%     if i > 1
%         box off
%         xlabel ''
%         ylabel ''
%         %set(gca, 'YTick', []);
%     end
%     %colorbar;
% end
% p.panel(5+numel(sp_), 'title', 'Shuffled', 'xlab', 'Fluctuation mode, {\iti}',...
%     'ylab', 'Number of cells', 'letter', char('c'+numel(sp_)+1), 'y_shift', 0.04);
% %subplot(1, numel(sp_)+1, numel(sp_)+1);
% mean_median_loadings_s = squeeze(mean(abs(res(1).median_loadings_s)));
% min_d = 30;
% ns = res(1).n_sizes;
% im_data = (mean_median_loadings_s(ns >= min_d,1:min_d));
% surf(1:min_d, ns(ns>=min_d), im_data, 'EdgeColor', 'none');
% view(2);
% %set(gca, 'YScale', 'log');
% set(gca, 'ColorScale', colorscale);
% caxis([0.03 0.25]);
% xlim([1 min_d]);
% ylim([min_d+10, 500]);
% %xlabel 'Fluctuation mode, i'
% %ylabel 'Number of cells'
% title('Shuffled', 'FontName', 'Helvetica', 'FontSize', 6, 'FontWeight', 'normal', 'Color', 'r');
% 
% set(gca, 'FontSize', 6);
% set(gca, 'FontName', 'Helvetica');
% set(gca, 'TickLength', [0.02 0.02]);
% %set(gca, 'YTick', [1 2 4 8 16 32 64].*min_d);
% box off
% %axis off
% xlabel ''; ylabel '';
% %set(gca, 'YTick', []);
% colorbar;
% %set(gcf, 'Units', 'inches');
% %set(gcf, 'Position', [8.5521    6.2292    8.3125    1.6146]);
% colormap parula;
% p.format;
%%
% p.panel([9 10], 'letter', 'h', 'y_shift', 0.04,...
%     'xlab', 'Number of cells', 'ylab', 'max_{\iti}|cos(PC_{\iti}, \Delta{\it\mu})|');
% MultiSessionVisualizer.plot_single_filtered(n_sizes, series, {'r', 'b'}, [1 2 3]);
% set(gca, 'XScale', 'log');
% set(gca, 'YScale', 'log');
% xlabel 'Number of cells'
% ylabel 'max_i|cos(PC_i, Dm)|'
% xlim([-Inf 500]);
% ylim([-Inf 1]);
% p.format;

%% 
save_file = sprintf('adjacent_metrics_no_decoder_%s.mat', ntype);
if ~exist(save_file, 'file')
    for i_s_i = 1:numel(selected_indices)
        my_ticker = tic;
        s_i = selected_indices(i_s_i);
        res__(i_s_i) = adjacent_decoder_noise_runner(s_i, true, ntype);
        toc(my_ticker)
        fprintf('done index %d\n', s_i);
    end
    save(save_file, 'res__');
end

load(save_file);

%TODO make the graphs onto the Pub, then copy all graphs over into the 
% new pres (based on the old pres). 
% Include a couple slides at the end about LTM extraction

cutoff = 100;
asymp_line = @(n,m) Utils.fitaline(n,m,cutoff);

res_lookup = @(res, code) arrayfun(@(x)cellfun(@(y)median(y{code(1),code(2:3)}), x.results_table),res,'UniformOutput',false);
n = {res__.n_sizes};
signal = res_lookup(res__, 'mse');
noise = res_lookup(res__, 'mce');
noise_shuf = res_lookup(res__, 'mue');

snr = cellfun(@rdivide, signal, noise, 'UniformOutput', false);
snr_shuf = cellfun(@rdivide, signal, noise_shuf, 'UniformOutput', false);

[signal_slope, signal_slope_conf] = cellfun(asymp_line, n, signal, 'UniformOutput', false);
[noise_slope, noise_slope_conf] = cellfun(asymp_line, n, noise, 'UniformOutput', false);
[noise_shuf_slope, noise_shuf_slope_conf] = cellfun(asymp_line, n, noise_shuf, 'UniformOutput', false);

[asymp_snr, asymp_snr_conf] = cellfun(@uncertain_divide, signal_slope, signal_slope_conf, noise_slope, noise_slope_conf);
[asymp_snr_shuf, asymp_snr_shuf_conf] = cellfun(@uncertain_divide, signal_slope, signal_slope_conf, noise_shuf_slope, noise_shuf_slope_conf);



dff2_lim = [0 24];

p.panel([11 12], 'y_shift', 0.04, 'letter', 'i', 'xlab', 'Number of cells', 'ylab', 'Signal ([event rate]^2)');
MultiSessionVisualizer.plot_single_filtered(n, {signal}, {'k'}, [true true true]);
ylim(dff2_lim);

p.panel([13 14], 'y_shift', 0.04, 'letter', 'j', 'xlab', 'Number of cells', 'ylab', 'Noise ([event rate]^2)');
MultiSessionVisualizer.plot_single_filtered(n, {noise_shuf, noise}, {'r', 'b'}, [true true true]);
%ylim(dff2_lim);
%text(100, 1/3, 'Real', 'Color', 'b');
%text(300, 0.3/3, 'Shuffled', 'Color', 'r');
p.format;

p.panel([15 16], 'y_shift', 0.04, 'letter', 'k', 'xlab', 'Number of cells', 'ylab', 'Signal / Noise');
%for i_ = 1:numel(asymp_snr)
%    shadedErrorBar(n{1}, repmat(asymp_snr(i_), size(n{1})), repmat(asymp_snr_conf(i_), size(n{1})), 'lineprops', ':b');
%    hold on;
%end
MultiSessionVisualizer.plot_single_filtered(n, {snr, snr_shuf}, {'b', 'r'}, [true true true]);
%text(250, 6.2, 'Shuffled', 'Color', 'r');
%text(400, 2, 'Real', 'Color', 'b');
hold on;
%l_ = refline(0, asymp_snr); l_.Color = 'b';


p.format;
%%
p.print('supplements_pdf', sprintf('on_%s.pdf', ntype));
%%
function [quotient, quotient_uncertainty] = uncertain_divide(x, xc, y, yc)
quotient = x./y;
quotient_uncertainty = abs(x./y)*sqrt((xc./x).^2 + (yc./y).^2);
end