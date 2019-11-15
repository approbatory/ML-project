function union_mask = compute_union_mask(ds)

cell_inds = find(ds.is_cell);
num_cells = length(cell_inds);

union_mask = false(size(ds.cells(1).im), 'logical');

for k = 1:num_cells
    cell_ind = cell_inds(k);
    union_mask = union_mask | ds.cells(cell_ind).mask;
end