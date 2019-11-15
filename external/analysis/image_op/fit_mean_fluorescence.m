function F_fit = fit_mean_fluorescence(F)
% Fit a sum of two exponentials to the mean flourescence F
%
% 2015 01 31 Tony Hyun Kim (Revised: Hakan Inan, 15-Jan-4)
%

num_frames = length(F);

time = 0:(num_frames-1);

f2 = fit(time',F','exp2'); % Fit sum of 2 exponentials
params = coeffvalues(f2);
F_fit = params(1)*exp(params(2)*time) + params(3)*exp(params(4)*time);