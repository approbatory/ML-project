function filter_for_size(ds, size_threshold)
% Automatically rejects candidate cells if the number of pixels occupied by
% the filter does not exceed a threshold.
num_filtered = 0;
for k = 1:ds.num_cells
    num_mask_pixels = sum(ds.cells(k).mask(:));
    if num_mask_pixels < size_threshold
        num_filtered = num_filtered + 1;
        ds.cells(k).label = 'not a cell';
    end
end
fprintf('%s: %d cells filtered by size criterion; %d remaining cells\n',...
    datestr(now), num_filtered, ds.num_cells - num_filtered);