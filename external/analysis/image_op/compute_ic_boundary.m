function [boundaries, ic_mask] = compute_ic_boundary(ic_filter, ic_filter_thresh)
% Computes the boundary of the IC filter. Boundaries are returned as a cell
%   array, in decreasing order of the boundary length

A = threshold_ic_filter(ic_filter, ic_filter_thresh);

ic_mask = (A>0);
boundaries = bwboundaries(ic_mask, 'noholes');

% Sort boundaries in order of decreasing length
boundaries_lengths = cellfun(@length, boundaries);
[~, sort_idx] = sort(boundaries_lengths, 'descend');
boundaries = boundaries(sort_idx);

% Arrange each boundary as [x y]
boundaries = cellfun(@fliplr, boundaries, 'UniformOutput', false);