function stats = compute_trace_stats(trace)
% Compute trace statistics, including:
%   - Basic histogram
%   - Mode
%   - Percentiles (95, 96, 97, 98, 99-th)
%

stats.num_bins = 500; % FIXME: Set dynamically
[n, x] = hist(trace, stats.num_bins);

nonzero_inds = n>0;
stats.hist_counts = n(nonzero_inds);
stats.hist_centers = x(nonzero_inds);

[~, max_ind] = max(n);
stats.mode = x(max_ind);

% Also compute percentiles
ps = [95 96 97 98 99]';
stats.percentiles = [ps prctile(trace,ps)];