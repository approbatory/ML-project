function [fitresult, gof] = fishinfo_fit(cell_nums, te_fish_mean_preshuf)
%CREATEFIT(CELL_NUMS,TE_FISH_MEAN_PRESHUF)
%  Create a fit.
%
%  Data for 'preshuf' fit:
%      X Input : cell_nums
%      Y Output: te_fish_mean_preshuf
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 03-Aug-2018 18:08:46


%% Fit: 'preshuf'.
[xData, yData] = prepareCurveData( cell_nums, te_fish_mean_preshuf );

% Set up fittype and options.
ft = fittype( 'a*x/(1+e*x)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.0901637920024458 0.0646918793241604];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
% 
% % Plot fit with data.
% figure( 'Name', 'preshuf' );
% h = plot( fitresult, xData, yData );
% legend( h, 'te_fish_mean_preshuf vs. cell_nums', 'preshuf', 'Location', 'NorthEast' );
% % Label axes
% xlabel cell_nums
% ylabel te_fish_mean_preshuf
% grid on

