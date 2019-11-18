%% Figure 1 | Noise correlations limit the place information encoded in the neural ensemble
% (a) Schematic of the experiment.
% (Top) Mouse runs back and forth on a linear track motivated by water rewards on each end. 
% (Bottom) Position traces during running behavior.
% We include in the analysis only the highly stereotyped running motions (black), and treat them as separate trials.

figure;
data_source = DecodeTensor.cons_filt(70, true);
load(data_source{1});
pos_raw = tracesEvents.position(:,1);
time = ((1:numel(pos_raw)) - 26000)/20 - 55;
pos_norm = (pos_raw - prctile(pos_raw, 5))./(prctile(pos_raw, 95) - prctile(pos_raw, 5))*120; %cm
pos_norm(pos_norm < 0) = 0;
pos_norm(pos_norm > 120) = 120;
plot(time, pos_norm, 'k');
opt = DecodeTensor.default_opt; opt.total_length = 120; opt.bin_width = 6;
[~, ~, tr_s, tr_e, ~, ~, ~] = DecodeTensor.new_sel(pos_raw, opt);
bad_y = pos_norm;
for tr_i = 1:numel(tr_s)
    bad_y(tr_s(tr_i):tr_e(tr_i)) = nan;
end
hold on;
h = plot(time, bad_y, 'r');
xlim([0 4030]/20);
ylabel 'Position (cm)'
xlabel 'Time (s)'
%legend(h, 'Excluded frames');
%legend boxoff
ylim([0 120]);
%line([139 147], [140 140], 'Color', 'r');
%text(149, 140, 'Excluded frames');
annotation('textbox', [(163.5/150-0.412) 0.75 0.3 0.3], 'String', 'Excluded frames',...
    'FitBoxToText', 'on', 'FontName', 'Helvetica', 'FontSize', 6, 'LineStyle', 'none');
annotation('line', [142 147]/150-0.295, [140 140]/150 + 0.038, 'Color', 'r')
%set(gca, 'XTick', []);
%set(gca, 'YTick', []);
figure_format('boxsize', [3.625 0.77], 'fig_factor', 1.45);
Utils.printto('figure1_pdf/demo', 'position_trajectory.pdf');

%% 
% (c) Example place decoding result from several trials, showing correct place bin (black),
% predicted place bin (blue), and prediction from shuffled data (red). Each place bin is 6 cm wide.

PanelGenerator.decode_demo(true);

%%
% (d) The rate at which the place decoder misidentifies a given place bin,
% as a function of its position along the track, shown for real data (blue)
% and shuffled data (red) from 497 neurons in one animal.
% Showing leftward and rightward runs separately. Errorbars represent SEM over random train-test divisions.
PanelGenerator.confusion(true);

%%
% (e) Inverse mean squared error of a neural decoder for place as a
% function of the number of cells included in the decoder.
% Each point is an aggregate of 80 random subsets of cells at fixed size.
% Shown for 3 mice. Errorbars represent SEM over random subsets and train-test divisions.
% (Inset) The equivalent curves as the square root of the mean squared error.
PanelGenerator.decoding_curves(true);

%%
% (f) Schematic of the equation for Fisher information as a function of ensemble size for an
% information-limited code (blue) and an equivalent decorrelated code (red).
% The parameter I0 represents linear slope at small ensemble size. 
% N is the ensemble size at which Fisher information is half compared to the decorrelated code.

figure;
x_vals = 1:500;
I_func = @(n) n./(1 + n./100);
plot(x_vals, I_func(x_vals), 'b'); hold on;
plot(x_vals, x_vals, 'r');
ylim([0 200]);
xlim([0 500]);
l_ = refline(0, 100); l_.Color = 'k';
set(gca, 'XTick', 0:100:500);
set(gca, 'XTickLabel', {'0', 'N', '2N', '3N', '4N', '5N'});
m_ = 1;
set(gca, 'YTick', [m_./(m_+1).*100, 100]);
set(gca, 'YTickLabel', {'I_0N/2', 'I_0N'});
%line([100 100], [0 500], 'LineStyle', ':', 'Color', 'k');
%line([0 500], [50 50], 'LineStyle', ':', 'Color', 'k');
text(200, 20, 'I_0n/(1+n/N)', 'Color', 'b');
text(200, 165, 'I_0n', 'Color', 'r');
text(350, 100, '99% I_0N', 'Color', 'g', 'VerticalAlignment', 'bottom');
l_ = refline(0, 99); l_.Color = 'g'; l_.LineStyle = ':';
xlabel 'Number of cells'
ylabel '1/MSE'
figure_format;
Utils.printto('figure1_pdf/demo', 'saturation_curve.pdf');

%%
% (g, h) Fitted parameters to the function in f, showing initial linear slope in g and ensemble size at half maximum in h.
% Values come from 107 sessions from 12 mice, shown in aggregate for each mouse.
% Colors denote mouse identity and errorbars are aggregated 95% confidence intervals.
% Initial linear slope shows no significant change (P > 0.05) but an increase in ensemble size at half maximum by an order of magnitude
% (P < 0.001) when comparing shuffled to unshuffled data using the Wilcoxon signed rank test.
% <plots produced under e>

%% Figure 2 | Noise correlations contain modes that remain aligned with the position tuning curve as ensemble size increases,
%% resulting in linearly rising noise variance in the coding direction
% (a) Two-dimensional projection of real (left) and shuffled (right) neural data using partial least squares regression on position and direction of motion.
% Color coded by place bin (bottom).
d = DecodeTensor.cons_filt(70); %DecodeTensor(6, 'IED');
[X, ks] = d.get_dataset;
X_shuf = shuffle(X, ks);
X_z = zscore(X);
X_shuf_z = zscore(X_shuf);
[XS, stats, origin] = Utils.pls_short(X_z, [ceil(ks/2), mod(ks,2)]);
XS_s = (X_shuf_z - mean(X_shuf_z)) * stats.W;
origin_s = -mean(X_shuf_z) * stats.W;

figure;
scatter(XS(:,1), XS(:,2), 8, Utils.colorcode(ceil(ks/2)), 'filled', 'MarkerFaceAlpha', 0.02); hold on;
scatter(origin(1), origin(2), 10, 'k');
xlabel PLS1; ylabel PLS2;
axis equal; axis tight;
xl_ = xlim;
yl_ = ylim;
set(gca, 'XTick', []);
set(gca, 'YTick', []);
figure_format('boxsize', [1 1.2], 'fontsize', 6);
print('-dpng', '-r1800', 'figure2_pdf/demo/PLS_adjacent.png');

figure;
scatter(XS_s(:,1), XS_s(:,2), 8, Utils.colorcode(ceil(ks/2)), 'filled', 'MarkerFaceAlpha', 0.02); hold on;
scatter(origin_s(1), origin_s(2), 10, 'k');
xlabel PLS1; ylabel PLS2;
axis equal; axis tight;
xlim(xl_); ylim(yl_);
set(gca, 'XTick', []);
set(gca, 'YTick', []);
figure_format('boxsize', [1 1.2], 'fontsize', 6);
print('-dpng', '-r1800', 'figure2_pdf/demo/PLS_adjacent_shuf.png');

%%
% (c) Absolute cosine of the angle between the principal noise correlation mode most aligned with the signal direction,
% as a function of ensemble size. Real data (blue) and shuffled data (red) from three mice. Showing the median over place bins.
% Shaded region represents 95% confidence interval over 20 random subsets of cells.
% (Top schematic) Illustration of the effect of a correlation mode remaining aligned with a specific direction as one adds a new dimension (blue).
% For isotropic noise the alignment will decay as the inverse square root of dimension (red).
% (Inset) The average slope of the curves on the log-log plot, computed beyond an ensemble size of 50 and representing a power-law exponent.
% Values for 107 sessions from 12 mice, aggregated within each mouse. Errorbars denote aggregated 95% confidence intervals.
% Color denotes mouse identity as in Figure 1g,h. Shuffling reduces the exponent (P < 0.001, Wilcoxon signed rank test).

%<depends on precomputed MedLoad_agg_190705-171806_0.mat>
PanelGenerator.medload;

%%
% (d) Projected noise variance and squared distance between mean responses (signal strength, black)
% for discrimination between adjacent place bins, as a function of ensemble size.
% Projecting real data (blue) and shuffled data (red) onto the direction of a correlation-insensitive classifier,
% and showing the median values over place bins. Values are normalized by the slope of the squared signal strength.
% Shown for three mice. Shaded regions are 95% confidence intervals over 20 random subsets of cells and train-test divisions for the classifier.

% 105 samples, not 107 because 2 sessions did not satisfy parameter
% confidences being < 0.5 * parameter values

%<depends on precomputed adjacent_agg_190725-194911_0.mat>
PanelGenerator.adjacent_decoders;

%%
% (e) The slope of the normalized noise variance from d for real and shuffled data, computed beyond ensemble sizes of 100.
% Values, errorbars, and color scheme obtained, aggregated, and presented as in c.
% Shuffling drastically reduces the slopes (P < 0.001, Wilcoxon signed rank test).
% (f) Linear slope coefficients I0, from Figure 1g, plotted against the inverse of the y-intercept value of the normalized noise variance from e, using shuffled data.
% Shown from 105 sessions from 12 mice, and color-coded by mouse identity. Errorbars are 95% confidence intervals from the fits.
% The variables are directly correlated (R2=0.47).
% (g) Asymptotic IMSE values I0N, computed from fitted parameters as in Figure 1g,h, plotted against the inverse of the slope of the normalized noise variance from e, using real data.
% Sessions, color scheme, and errorbars as in f. The variables are directly correlated (R2=0.43).

%< plots generated under d>

%% Figure 3 | Noise variance in the spatial encoding direction increases linearly with ensemble size.
signal_and_noise_analysis;

% (a) 