%comparing distances between cells to correlation coefficients
function neuropil_contamination_argument(neural_data_type)
d = DecodeTensor(4, neural_data_type);
[X, ks] = d.get_dataset();
[n_cells, n_place_bins, n_trials] = size(d.data_tensor);
%% presumably, neuropil noise will not be coding dependent,
% and so the contamination noise will not depend on signal
% so it can be calculated by aggregating over all bins
bin_list = unique(ks);
bin_mean_response = zeros(numel(bin_list), n_cells);
for i_b = 1:numel(bin_list)
    b = bin_list(i_b);
    bin_mean_response(i_b, :) = mean(X(ks == b,:));
end

copied_means = bin_mean_response(ks, :);
X_fluctuations = X - copied_means;
%% noise correlation matrix
rho = corr(X_fluctuations);
rho_signal = corr(bin_mean_response);
rho_total = corr(X);

right_mean_response = bin_mean_response(1:2:end,:);
left_mean_response = bin_mean_response(2:2:end,:);

[right_com, right_ord] = find_neuron_com(right_mean_response);
[left_com, left_ord] = find_neuron_com(left_mean_response);

%% distance matrix
pack = load(d.source_path);
cell_loc = pack.tracesEvents.cellAnatomicLocat;
%% calculate distances
cell_dist = sqrt((cell_loc(:,1) - cell_loc(:,1)').^2 + (cell_loc(:,2) - cell_loc(:,2)').^2);
r_vec = rho(:); d_vec = cell_dist(:); rs_vec = rho_signal(:);
% %% plot scatter of correlation coeffs and distances
% r_vec = rho(:); d_vec = cell_dist(:);
% figure;
% basic_filt = (r_vec < 1) & (d_vec > 0);
% scatter(r_vec(basic_filt), d_vec(basic_filt), 9, 'filled', 'MarkerFaceAlpha', 0.02);
% 
% [linfit_tot, goodness_tot] = fit(r_vec(basic_filt), d_vec(basic_filt), 'poly1');
% text(-0.75, 125, sprintf('adj. R^2 = %.2f', goodness_tot.adjrsquare));
% hold on;
% plot(linfit_tot);
% 
% %figure;
% r_cut = 0.1; d_cut = 30;
% line([r_cut r_cut], [0 d_cut], 'Color', 'black');
% line([r_cut 1], [d_cut d_cut], 'Color', 'black');
% rest_filt = (r_vec > r_cut) & (d_vec < d_cut) & basic_filt;
% %scatter(r_vec(rest_filt), d_vec(rest_filt), 9, 'filled', 'MarkerFaceAlpha', 0.1);
% %xlabel 'Noise correlation \rho'
% %ylabel 'Cell distance (pix)'
% [linfit, goodness] = fit(r_vec(rest_filt), d_vec(rest_filt), 'poly1');
% text(0.4, 40, sprintf('adj. R^2 = %.4f', goodness.adjrsquare));
% hold on;
% l_ = plot(linfit); 
% fit_plot_filt = l_.XData > r_cut;
% l_.XData = l_.XData(fit_plot_filt);
% l_.YData = l_.YData(fit_plot_filt);
% legend off;
% xlabel 'Noise correlation \rho'
% ylabel 'Cell distance (pix)'
%%
%sc_f(r_vec, d_vec, basic_filt & ~rest_filt, rest_filt,...
%    'Noise corr.', 'Cell distance (pix)', -0.75, 125, 0.4, 40);
figure;
subplot(1,3,1);
basic_filt = (r_vec < 1) & (d_vec > 0);
r_cut = 0.1; d_cut = 25;
rest_filt = (r_vec > r_cut) & (d_vec < d_cut) & basic_filt;
sc_f(d_vec, r_vec, basic_filt & ~rest_filt, rest_filt,...
    'Cell distance (pix)', 'Noise corr.', 140, 0.1, 25, 0.8);
ylim([-1 1]); xlim([0 300]);
% signal corr vs. noise corr
subplot(1,3,2);
sc_f(rs_vec, r_vec, basic_filt & ~rest_filt, rest_filt,...
    'Signal corr.', 'Noise corr.', 0, -0.1, 0, 0.8);

% signal corr vs. distance
subplot(1,3,3);
sc_f(rs_vec, d_vec, basic_filt & ~rest_filt, rest_filt,...
    'Signal corr.', 'Cell distance (pix)', -0.75, 125, -0.9, 10);

set(gcf, 'Position', [203 521 1409 367]);
end

function [com, ord] = find_neuron_com(mean_response)
n = size(mean_response, 1);
com = sum(mean_response./sum(mean_response) .* (1:n)');
[~, ord] = sort(com);
end

function sc_f(a, b, filt1, filt2, lab1, lab2, l1x, l1y, l2x, l2y)
hold on;
scatter(a(filt1), b(filt1), 9, 'filled', 'MarkerFaceAlpha', 0.02);
scatter(a(filt2), b(filt2), 9, 'filled', 'red', 'MarkerFaceAlpha', 0.02);
[linfit1, goodness1] = fit(a(filt1), b(filt1), 'poly1');
[linfit2, goodness2] = fit(a(filt2), b(filt2), 'poly1');
l_ = plot(linfit1); l_.Color = 'b';
l_ = plot(linfit2); l_.Color = 'r';
text(l1x, l1y, sprintf('adj. R^2 = %.2f', goodness1.adjrsquare));
text(l2x, l2y, sprintf('adj. R^2 = %.2f', goodness2.adjrsquare));
xlabel(lab1); ylabel(lab2);
legend off;
end