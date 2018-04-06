num_filtered = 0;
for k = 1:ds.num_cells
    num_mask_pixels = sum(ds.cells(k).mask(:));
    if num_mask_pixels < 50
        num_filtered = num_filtered + 1;
        ds.cells(k).label = 'not a cell';
    end
end
