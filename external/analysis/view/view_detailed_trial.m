function view_detailed_trial(ds, cell_idx, trial_idx)
% Examine the trace of a cell alongside the behavior video for that trial.
% The DaySummary must have the behavior video loaded.

prev_trial = 0;

if ~ds.is_behavior_loaded
    fprintf('  Behavior video has not been loaded into this DaySummary!\n');
else
    while (1)
        if is_valid_trial(trial_idx)
            if (trial_idx ~= prev_trial) % Don't need to redraw if cell_idx didn't change
                draw_cell_trial(cell_idx, trial_idx);
                prev_trial = trial_idx;
            end
        else
            fprintf('  Invalid trial index (%d) for this DaySummary. Exiting!\n', trial_idx);
            break;
        end
        
        % Ask user for command
        prompt = sprintf('Showing Cell %d / Trial %d >> ', cell_idx, trial_idx);
        resp = strtrim(input(prompt, 's'));
        
        val = str2double(resp);
        if ~isnan(val) % Is a number
            trial_idx = val;
        else % Not a number
            resp = lower(resp);
            switch (resp)
                case {'n', ''} % Next trial
                    if (trial_idx < ds.num_trials)
                        trial_idx = trial_idx + 1;
                    else
                        fprintf('  Already at final trial (%d)!\n', ds.num_trials);
                    end
                    
                case 'p' % Previous trial
                    if (1 < trial_idx)
                        trial_idx = trial_idx - 1;
                    else
                        fprintf('  Already at first trial!\n');
                    end
                    
                case 'q' % Exit
                    break;
                    
                otherwise
                    fprintf('  Could not parse "%s"\n', resp);
            end
        end
    end
end

    function valid_trial = is_valid_trial(trial_idx)
        valid_trial = (1 <= trial_idx) && (trial_idx <= ds.num_trials);
    end % is_valid_trial

    function draw_cell_trial(cell_idx, trial_idx)
        clf;
        f = gcf;
        
        Mb = ds.get_behavior_trial(trial_idx); % Behavior movie
        trace = ds.trials(trial_idx).traces(cell_idx,:);
        num_frames_in_trial = length(trace);
        
        % Image of cell
        %------------------------------------------------------------
        subplot(3,4,[1 2]);
        imagesc(ds.cells(cell_idx).im);
        axis image;
        title(sprintf('Cell %d (%s)', cell_idx, ds.cells(cell_idx).label));
        colormap jet; freezeColors;
        
        % Zoom into the selected trial in the full raster
        %------------------------------------------------------------
        trial_window = 10;
        trial_window_start = max(1, trial_idx-trial_window);
        trial_window_end = min(ds.num_trials, trial_idx+trial_window);
        
        subplot(3,4,[5 6 9 10]);
        ds.plot_cell_raster(cell_idx, 'draw_correct');
        scale = get(gca, 'CLim');
        ylim([trial_window_start-0.5 trial_window_end+0.5]);
        set(gca, 'YTick', trial_window_start:trial_window_end);
        x_ends = get(gca, 'XLim');
        line(x_ends, (trial_idx-0.5)*[1 1], 'Color', 'w', 'LineWidth', 2);
        line(x_ends, (trial_idx+0.5)*[1 1], 'Color', 'w', 'LineWidth', 2);
        
        title('All trials');
        colormap jet; freezeColors;
        
        % Show trace
        %------------------------------------------------------------
        trial_frames = double(ds.trial_indices(trial_idx, :)); % [start open-gate close-gate end]
        trial_frames = trial_frames - trial_frames(1) + 1;
        
        align_index = 3; % Align to closing of gate
        align_str = 'Frames relative to gate close';
        alignment_frame = trial_frames(align_index);
        
        trial_frames = trial_frames - alignment_frame; % Time relative to gate close
        
        subplot(3,4,[3 4]);
        plot(trial_frames(1):trial_frames(4), trace, 'LineWidth', 2, 'HitTest', 'off');
        xlim(trial_frames([1 4]));
        ylim(scale);
        grid on;
        title(sprintf('Trial %d', trial_idx));
        xlabel(align_str);
        ylabel('Signal [a.u.]');
        hold on;
        
        plot(trial_frames(2)*[1 1], scale, 'r--', 'HitTest', 'off'); % Open-gate
        plot(trial_frames(3)*[1 1], scale, 'r--', 'HitTest', 'off'); % Close-gate
        
%         if ds.is_tracking_loaded
%             trial = ds.trials(trial_idx);
%             movement_onset = trial.movement_onset_frame / num_trials_in_frame;
%             turn_onset = trial.turn_onset_frame / num_trials_in_frame;
%             plot(movement_onset*[1 1], scale, 'k--', 'HitTest', 'off');
%             plot(turn_onset*[1 1], scale, 'k--', 'HitTest', 'off');
%         end
        
        % Markers for indicating current frame
        t = plot(trial_frames(1)*[1 1], scale, 'ButtonDownFcn', @start_drag); % Vertical bar -- can be dragged
        d = plot(trial_frames(1), trace(1), 'o',... % Dot
                 'MarkerFaceColor', 'b', 'MarkerSize', 8,...
                 'HitTest', 'off');
        trace_axis = gca;
        set(trace_axis, 'ButtonDownFcn', @update_frame);
        set(f, 'WindowButtonUpFcn', @end_drag);
        
        % Show behavior movie
        %------------------------------------------------------------
        subplot(3,4,[7 8 11 12]);
        hb = imagesc(Mb(:,:,1));
        set(gca, 'XTick', []);
        set(gca, 'YTick', []);
        axis image;
        colormap gray;
        title(sprintf('Frame 1 of %d', num_frames_in_trial));
        
        % If tracking data loaded, overlay the positional information
        if (ds.is_tracking_loaded)
            hold on;
            centroids = ds.trials(trial_idx).centroids;
            plot(centroids(:,1), centroids(:,2), '.-');
            htrack = plot(centroids(1,1), centroids(1,2), 'ro');
        end
        
        function start_drag(~, ~)
            set(f, 'WindowButtonMotionFcn', @update_frame);
        end
        
        function end_drag(~, ~)
            set(f, 'WindowButtonMotionFcn', '');
        end
        
        function update_frame(~, ~)
            cp = get(trace_axis, 'CurrentPoint');
            sel_frame = round(cp(1)); % X point of click
            
            sel_frame = sel_frame + alignment_frame;
            sel_frame = max(sel_frame, 1);
            sel_frame = min(sel_frame, num_frames_in_trial);
            
            % Update visuals
            set(t, 'XData', (sel_frame-alignment_frame)*[1 1]);
            set(d, 'XData', (sel_frame-alignment_frame), 'YData', trace(sel_frame));
            set(hb, 'CData', Mb(:,:,sel_frame));
            
            subplot(3,4,[7 8 11 12]);
            title(sprintf('Frame %d of %d', sel_frame, num_frames_in_trial));
            
            if (ds.is_tracking_loaded)
                centroid = ds.trials(trial_idx).centroids(sel_frame, :);
                set(htrack, 'XData', centroid(1), 'YData', centroid(2));
            end
        end % update_frame
    end % draw_cell_trial

end % view_detailed_trial