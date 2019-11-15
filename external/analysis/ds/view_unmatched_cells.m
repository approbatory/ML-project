function unmatched_ids = view_unmatched_cells(ds1, m_1to2)

[unmatched_ids, matched_ids] = find_unmatched_cells(ds1, m_1to2);
num_unmatched_ids = length(unmatched_ids);

ds1.plot_cell_map({matched_ids, 'g'; unmatched_ids, 'y'});
title(sprintf('Found %d unmatched cells (yellow) out of %d cells',...
    num_unmatched_ids, ds1.num_classified_cells));