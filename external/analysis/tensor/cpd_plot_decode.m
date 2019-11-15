function [ln, pt] = cpd_plot_decode(models, depvar, best_n)
% PLOT_DECODE, plots decoding accuracy as a function of model rank for a list of models
%
%     plot_decode(models, 'start')
%     plot_decode(models, 'end')
%     plot_decode(models, 'correct')

if nargin == 2
	best_n = 1;
end

% collect rank and error of each model

[num_starts, max_rank] = size(models);

accuracy = nan(num_starts, max_rank);
ranks = nan(num_starts, max_rank);
for r = 1:max_rank
	for s = 1:best_n
	    ranks(s, r) = r;
	    accuracy(s, r) = models(s, r).decode.(depvar);
	end
end
ranks = ranks(:);
accuracy = accuracy(:);

unique_rank = unique(ranks(~isnan(ranks)));
n_unique = length(unique_rank);

% average err for each rank
avg_acc = zeros(n_unique,1);
for i = 1:n_unique
    avg_acc(i) = mean(accuracy(ranks == unique_rank(i)));
end

% make figure
hold on
ln = plot(unique_rank, avg_acc, '-r', 'linewidth', 2);
pt = plot(ranks, accuracy, '.k', 'markersize', 10);
title(depvar);
