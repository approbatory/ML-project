function [secondary_ids1, secondary_ids2] = view_secondary_cells(ds1, ds2, m_1to2, m_2to1)

% Transfer (positive) labels from ds2 --> ds1
%------------------------------------------------------------
is_cell1 = find(ds1.is_cell);
secondary_ids1 = find_secondary_cells(ds1, ds2, m_1to2);
is_not_cell1 = setdiff(find(~ds1.is_cell), secondary_ids1);
num_secondary1 = length(secondary_ids1);

subplot(121);
ds1.plot_cell_map({is_cell1, 'g'; secondary_ids1, 'c'; is_not_cell1, 'r'});
title(sprintf('Dataset 1: %d secondary cells labels (cyan)', num_secondary1));

% Transfer (positive) labels from ds1 --> ds2
%------------------------------------------------------------
is_cell2 = find(ds2.is_cell);
secondary_ids2 = find_secondary_cells(ds2, ds1, m_2to1);
is_not_cell2 = setdiff(find(~ds2.is_cell), secondary_ids2);
num_transfer_to_2 = length(secondary_ids2);

subplot(122);
ds2.plot_cell_map({is_cell2, 'g'; secondary_ids2, 'c'; is_not_cell2, 'r'});
title(sprintf('Dataset 2: %d transferred labels (cyan)', num_transfer_to_2));