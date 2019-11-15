function [unmatched_ids, matched_ids] = find_unmatched_cells(ds1, m_1to2)
% Find classified cells in provided DaySummary ('ds1') with no match in
% the corresponding matching ('m_1to2').
%
% NOTE: The result is dependent on _how_ the matching ('m_1to2') is
% computed. There are two possibilities:
%
%   (1) If the matching is computed in a way that omits non-cells (default
%       behavior of `match_masks`) then cells in ds1 that overlap with
%       non-cells in ds2 are considered non-matched.
%
%   (2) If the matching is computed including the non-cells, cells in ds1
%       are considered unmatched if they do not match to any source (cell
%       or non-cell) in ds2. This is a more strict definition of
%       "unmatched" than the first case, and will yield a smaller set.
%

is_cell = find(ds1.is_cell);

m_1to2 = m_1to2(is_cell);
is_not_matched = cellfun(@isempty, m_1to2);

unmatched_ids = is_cell(is_not_matched);
matched_ids = setdiff(is_cell, unmatched_ids);