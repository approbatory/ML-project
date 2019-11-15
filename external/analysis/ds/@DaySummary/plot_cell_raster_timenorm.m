function raster = plot_cell_raster_timenorm(obj, cell_idx, varargin)
    % Optional argument 'draw_correct' will place a box at the end
    % of each trial indicating correct (green) or incorrect (red)
    %
    % Additional optional arguments allow for filtering of trials,
    % e.g. "plot_cell_raster(cell_idx, 'start', 'east')"

    display_trial = ones(obj.num_trials, 1);
    draw_correct = 0;
    if ~isempty(varargin)
        for k = 1:length(varargin)
            vararg = varargin{k};
            if ischar(vararg)
                switch lower(vararg)
                    case 'draw_correct'
                        draw_correct = 1;
                end
            end
        end
        % Trial filtering arguments
        display_trial = obj.filter_trials(varargin{:});
    end

    resample_grid = linspace(0, 1, 1000);
    num_filtered_trials = sum(display_trial);
    raster = zeros(num_filtered_trials, length(resample_grid));
    correctness = zeros(num_filtered_trials);
    counter = 0;
    for k = 1:obj.num_trials
        if display_trial(k)
            counter = counter+1;
            line = obj.get_trace(cell_idx, k);
            raster(counter,:) = interp1(linspace(0,1,length(line)),...
                              line,...
                              resample_grid,...
                              'pchip');
            correctness(counter) = obj.trials(k).correct;
        end
    end

    imagesc(resample_grid, 1:num_filtered_trials, raster);
    colormap jet;
    xlabel('Trial phase [a.u.]');
    xlim([0 1]);
    ylabel('Trial index');
    set(gca, 'TickLength', [0 0]);
    
    if (draw_correct)
        corr_width = 0.025;
        for k = 1:num_filtered_trials
            if correctness(k)
                corr_color = 'g';
            else
                corr_color = 'r';
            end
            rectangle('Position', [1 k-0.5 corr_width 1],...
                      'FaceColor', corr_color);
        end
        xlim([0 1+corr_width]);
    end
end