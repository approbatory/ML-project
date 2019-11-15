function view_detailed_raster(ds, cell_idx)
% Plot a detailed "raster summary" of a single cell in DaySummary
%
% Usage:
%     m3d15 = DaySummary(sources.maze, 'rec001', 'reconst', 'excludeprobe');
%       14-Aug-2015 13:58:23: Loaded data from rec001/rec_150813-115003.mat
%       14-Aug-2015 13:58:23: Loaded classification from rec001/class_150813-115320.txt
%     browse_rasters(m3d15);
%

while (1)
    draw_rasters(); % Global normalized scaling for all trials
    
    % Ask user for command
    prompt = sprintf('Cell %d raster >> ', cell_idx);
    resp = strtrim(input(prompt, 's'));
    
    val = str2double(resp);
    if (~isnan(val)) % Is a number (trial index)
        if ds.is_behavior_loaded
            if (val == 0) % Trial "0" shows all trials superimposed in space
                if ds.is_tracking_loaded
                    view_trace_in_space(ds, cell_idx);
                    fprintf('  Showing cell %d trace on space. Press any key to continue...\n', cell_idx);
                    pause;
                else
                    fprintf('  Tracking data not loaded into DaySummary!\n');
                end
            else
                view_detailed_trial(ds, cell_idx, val);
            end
        else
            fprintf('  Behavior video not loaded into DaySummary!\n');
        end
    else
        resp = lower(resp);
        switch (resp)
            case 't' % "Take" screenshot
                savename = sprintf('raster_%d.png', cell_idx);
                print('-dpng', savename);
                fprintf('  Saved raster as "%s"!\n', savename);
            
            case {'q', ''} % Exit
                break;
                
            otherwise
                fprintf('  Could not parse "%s"\n', resp);
        end
    end
end

    % Helper functions
    %------------------------------------------------------------
    function draw_rasters()
        clf;
        
        % Image of cell
        subplot(3,4,[1 2]);
        imagesc(ds.cells(cell_idx).im);
        axis image;
        title(sprintf('Cell %d (%s)', cell_idx, ds.cells(cell_idx).label));
        colormap jet; freezeColors;
        
        % Raster of all trials, with correctness
        subplot(3,4,[5 6 9 10]);
        ds.plot_cell_raster(cell_idx, 'draw_correct');
        raster_scale = get(gca, 'CLim'); % Scale that applies to all trials
        title('All trials');
        colormap jet; freezeColors;
        
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
    end % draw_rasters
end % view_cell_rasters