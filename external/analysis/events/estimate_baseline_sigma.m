function [baseline, sigma, stats] = estimate_baseline_sigma(trace)

stats = compute_trace_stats(trace);

baseline = stats.mode;
tr_lower = trace(trace <= baseline);
sigma = std(tr_lower - baseline);

% Convert standard deviation of the half-normal distribution 
% to that of the normal distribution
sigma = sigma / (1 - 2/pi);
% thresh = stats.mode + 3*sigma;
