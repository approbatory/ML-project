function draw_standard_subrasters(ds, cell_idx, raster_scale)
% Used by 'view_raster_touch'

% Divide rasters by correctness
subplot(3,4,3);
ds.plot_cell_raster(cell_idx, 'correct');
set(gca, 'CLim', raster_scale);
title('Correct');
subplot(3,4,4);
ds.plot_cell_raster(cell_idx, 'incorrect');
set(gca, 'CLim', raster_scale);
title('Incorrect');

% Divide rasters by start location
subplot(3,4,7);
ds.plot_cell_raster(cell_idx, 'start', 'west', 'draw_correct');
set(gca, 'CLim', raster_scale);
title('West start'); 
subplot(3,4,8);
ds.plot_cell_raster(cell_idx, 'start', 'east', 'draw_correct');
set(gca, 'CLim', raster_scale);
title('East start');

% Divide rasters by end location
subplot(3,4,11);
ds.plot_cell_raster(cell_idx, 'end', 'south', 'draw_correct');
set(gca, 'CLim', raster_scale);
title('South end');
subplot(3,4,12);
ds.plot_cell_raster(cell_idx, 'end', 'north', 'draw_correct');
set(gca, 'CLim', raster_scale);
title('North end');

end % draw_standard_subrasters