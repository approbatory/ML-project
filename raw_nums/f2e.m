function f2e
%% Figure 2e
% Depends on file: decoding_curves_fits.mat
% Raw numbers for:
% - LineFitX, LineEx1X, LineEx2X, LineEx3X
% - LineRealFitEx1Y, LineRealEx1Y, ErrorBarsRealEx1Y
% - LineShuffledFitEx1Y, LineShuffledEx1Y, ErrorBarsShuffledEx1Y
% - LineRealFitEx2Y, LineRealEx2Y, ErrorBarsRealEx2Y
% - LineShuffledFitEx2Y, LineShuffledEx2Y, ErrorBarsShuffledEx2Y
% - LineRealFitEx3Y, LineRealEx3Y, ErrorBarsRealEx3Y
% - LineShuffledFitEx3Y, LineShuffledEx3Y, ErrorBarsShuffledEx3Y
rng(0);

data = rms_decoding_curves;
make_xlsx(data, 'f2e');
end

%% helper functions
function data = rms_decoding_curves
savedir = 'events_figs/f2_events/decoding_curves';
ap = @(x) fullfile(savedir, x);

load 'decoding_curves_fits.mat' ...
    sess mouse_names n_sizes imse imse_s I0_fit I0_fit_s N_fit N_fit_s ...
    series_gof series_exp_gof;

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


fname = ap('decoding_curve_fit.pdf');
data = my_aux_decoding_curves(fname, sess, mouse_names, n_sizes, imse, imse_s,...
    I0_fit, I0_fit_s, N_fit, N_fit_s, 'b', 'r',...
    [0 0.16], [2 50], [0 Inf], true);
end

function data = my_aux_decoding_curves(fname, sess, mouse_names, n_sizes, imse, imse_alt,...
    I0_fit, I0_fit_alt, N_fit, N_fit_alt, color, color_alt,...
    ylim_imse, ylim_rms, ylim_multi, non_inset)
%For two series of IMSE curves (i.e. unshuffled-shuffled, or
%unshuffled-diagonal), plot selected sessions as IMSE, RMS, along with the fit curve and
%plot all sessions by mouse.
%Inputs:
%   fname: base filename with directory,
%       'figure1_pdf/decoding_curves/fit_decoding_curves.pdf'
%   sess: a char cell array of session IDs
%   mouse_names: char cell array of mouse names corresponding
%       to session IDs in sess
%   n_sizes: x-axis values for each session in sess, as cell
%   imse: y-axis values for each session in sess, as cell
%   imse_alt: y-axis values for each session in sess, as cell,
%       for the other condition
%   I0_fit: I0 fit value for each session in sess, as numeric array
%   I0_fit_alt: I0 fit value for each session in sess, as
%       numeric array, for the other condition
%   N_fit: N fit value for each session in sess, as numeric
%       array
%   N_fit_alt: N fit value for each session in sess, as numeric
%       array, for the other condition
%   color: The color of the curves, e.g. 'b', 'r', 'm' (should
%       be 'b' for unshuffled)
%   color_alt: The color of the curves, for the other condition
%       e.g. 'b', 'r', 'm' (should be 'r' for shuffled, 'm' for diagonal)
%   ylim_imse: The y axis limits for the IMSE curves
%   ylim_rms: The y axis limits for the RMS curves
%   ylim_multi: The common y axis limits for the complete set
%       of curves. If given as [y1 y2] then y1 is ignored and
%       assumed to be 0. Otherwise sets only the upper y limit.
%       set to Inf for each mouse to have its own max
%figure('FileName', fname);
show_mice = {'Mouse2022', 'Mouse2024', 'Mouse2028'};

sp_ = SessManager.special_sessions_usable_index(show_mice);

%{
data_shuf = my_plot_decoding_curve(sess, sp_, n_sizes, imse_alt, I0_fit_alt, N_fit_alt, color_alt);
hold on;
data_real = my_plot_decoding_curve(sess, sp_, n_sizes, imse, I0_fit, N_fit, color);
xlabel 'Number of cells'
ylabel '1/MSE (cm^{-2})'

ylim(ylim_imse);
legend off
figure_format([2 2.5]);

data.LineFitX = data_real(1).LineFitX;
%}

[pref, fn, ext] = fileparts(fname);
figure('FileName', fullfile(pref, ['inset' fn ext]));
data_shuf = my_plot_decoding_curve(sess, sp_, n_sizes, imse_alt, I0_fit_alt, N_fit_alt, color_alt, true);
hold on;
data_real = my_plot_decoding_curve(sess, sp_, n_sizes, imse, I0_fit, N_fit, color, true);
l_ = refline(0, 5); l_.Color = 'k'; l_.LineStyle = ':';
xlabel 'Number of cells'
ylabel 'RMS Error (cm)'

ylim(ylim_rms);
set(gca, 'YScale', 'log');
set(gca, 'YTick', [1 2 5 10 20 50]);
legend off
if exist('non_inset', 'var') && non_inset
    figure_format([2 2.5]);
else
    figure_format('boxsize', [0.6 0.8], 'fontsize', 5);
end

data.LineFitX = data_real(1).LineFitX;
perm = [2 3 1];
for jx = 1:3
    ix = perm(jx);
    data.("LineEx"+jx+"X") = data_real(ix).LineX;
    
    data.("LineRealFitEx"+jx+"Y") = data_real(ix).LineFitY;
    data.("LineRealEx"+jx+"Y") = data_real(ix).LineY;
    data.("ErrorBarsRealEx"+jx+"Y") = data_real(ix).ErrorBarsY;
    
    data.("LineShuffledFitEx"+jx+"Y") = data_shuf(ix).LineFitY;
    data.("LineShuffledEx"+jx+"Y") = data_shuf(ix).LineY;
    data.("ErrorBarsShuffledEx"+jx+"Y") = data_shuf(ix).ErrorBarsY;
end

%MultiSessionVisualizer.plot_series(n_sizes, {imse_alt, imse}, {color_alt, color}, mouse_names, ylim_multi);
%xlabel 'Number of cells'
%ylabel '1/MSE (cm^{-2})'
%multi_figure_format;
end

function data = my_plot_decoding_curve(sess, sp_, n_sizes, imse_s, I0_fit_s, N_fit_s, color, isrms)
%plot a subset of the sessions as decoding curves + curve fits
%Inputs:
%   sess: cell array of strings of the session codes
%   sp_: subset of the indices in sess to display
%   n_sizes: x-axis values for each session in sess, as cell
%   imse_s: y-axis values for each session in sess, as cell
%   I0_fit_s: I0 fit value for each session in sess, as numeric array
%   N_fit_s: N fit value for each session in sess, as numeric
%       array
%   color: The color of the curves, e.g. 'b', 'r', 'm'
%   isrms (optional, default false): plot as RMS instead of
%      IMSE, this applies @(x)x.^(-1/2) to the values in imse_s
if ~exist('isrms', 'var')
    isrms = false;
end
if ~isrms
    pre = @(x)x;
    cap = 1;
else
    pre = @(x)x.^(-1/2);
    cap = 0.4;
end

ix = 0;
for j = 1:numel(sess)
    if ismember(j, sp_)
        ix = ix + 1;
        n = n_sizes{j};
        %i = imse{i};
        i_s = imse_s{j};
        %plot(series_fits{2}{j}, 'r');
        n_f = 1:500;
        plot(n_f, pre(I0_fit_s(j).*n_f./(1 + n_f./N_fit_s(j))), color);
        hold on;
        
        data(ix).LineFitX = n_f;
        data(ix).LineFitY = pre(I0_fit_s(j).*n_f./(1 + n_f./N_fit_s(j)));
        
        errorbar(n, mean(pre(i_s)), std(pre(i_s))./sqrt(size(i_s,1)),...
            color, 'Capsize', cap, 'LineStyle', 'none');
        
        data(ix).LineX = n;
        data(ix).LineY = mean(pre(i_s));
        data(ix).ErrorBarsY = std(pre(i_s))./sqrt(size(i_s,1));
    end
end
end