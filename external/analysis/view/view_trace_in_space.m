function view_trace_in_space(ds, cell_idx)
% Displays the trace of 'cell_idx' over the animal's trajectory in the
% behavior arena

[m, M] = get_trace_bounds(ds, cell_idx);
raster_scale = [m M];

clf;

subplot(2,4,1);
ds.plot_cell_raster(cell_idx, 'start', 'west', 'draw_correct');
set(gca, 'CLim', raster_scale);
title('West start');

subplot(2,4,8);
ds.plot_cell_raster(cell_idx, 'start', 'east', 'draw_correct');
set(gca, 'CLim', raster_scale);
title('East start');

subplot(2,4,5);
ds.plot_cell_raster(cell_idx, 'end', 'south', 'draw_correct');
set(gca, 'CLim', raster_scale);
title('South end');

subplot(2,4,4);
ds.plot_cell_raster(cell_idx, 'end', 'north', 'draw_correct');
set(gca, 'CLim', raster_scale);
title('North end'); 

subplot(2,4,[2 3 6 7]);
bg_image = ds.behavior_ref_img;
bg_image = cat(3, bg_image, bg_image, bg_image);
imagesc(bg_image);
axis image;
hold on;
title(sprintf('Cell %d (%s)', cell_idx, ds.cells(cell_idx).label));

for trial_idx = 1:ds.num_trials
    centroids = ds.trials(trial_idx).centroids;
    trace = ds.get_trace(cell_idx, trial_idx);
    trace = 255*(trace-m)/(M-m); % Scale trace to image colormap
    cline(centroids(:,1), centroids(:,2), trace, trace);
end
hold off;

end % view_trace_in_space