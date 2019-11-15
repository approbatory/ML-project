function secondary_ids = find_secondary_cells(ds1, ds2, m_1to2)
% Find sources in 'ds1' that are not cells, but their corresponding source
% in 'ds2' (under 'm_1to2') are classified cells.

is_not_cell = find(~ds1.is_cell);
is_cell_in_2 = false(size(is_not_cell));

for k = 1:length(is_not_cell)
    cell_idx1 = is_not_cell(k);
    match = m_1to2{cell_idx1};
    if ~isempty(match)
        cell_idx2 = match(1);
        is_cell_in_2(k) = ds2.is_cell(cell_idx2);
    end
end

secondary_ids = is_not_cell(is_cell_in_2);