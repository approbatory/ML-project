function view_detailed_trial_touch(ds, cell_idx, trial_idx, varargin)

subraster_type = '';
h_fig = [];
for i = 1:length(varargin)
    vararg = varargin{i};
    if ischar(vararg)
        switch lower(vararg)
            case 'fig'
                h_fig = varargin{i+1};
            case 'type'
                % The subraster type that should be requested when
                % "returning" to view_raster_touch
                subraster_type = varargin{i+1};
        end
    end
end

if isempty(h_fig)
    h_fig = figure;
else
    figure(h_fig);
    clf;
end

% Trial-specific data
%------------------------------------------------------------
[trial_frames, alignment_frame] = get_frames_relative_to_close(ds, trial_idx);
Mb = ds.get_behavior_trial(trial_idx); % Behavior movie
trace = ds.trials(trial_idx).traces(cell_idx,:);
num_frames_in_trial = length(trace);

% Image of cell
%------------------------------------------------------------
subplot(3,2,1);
imagesc(ds.cells(cell_idx).im);
axis image;
title(sprintf('Cell %d (%s)', cell_idx, ds.cells(cell_idx).label));
colormap jet; freezeColors;

% Zoom into the selected trial in the full raster
%------------------------------------------------------------
trial_window = 10;
trial_window_start = max(1, trial_idx-trial_window);
trial_window_end = min(ds.num_trials, trial_idx+trial_window);
gui.raster = subplot(3,2,[3 5]);
ds.plot_cell_raster(cell_idx, 'draw_correct', 'draw_events');
trace_scale = get(gca, 'CLim');
ylim([trial_window_start-0.5 trial_window_end+0.5]);
x_ends = get(gca, 'XLim');
line(x_ends, (trial_idx-0.5)*[1 1], 'Color', 'w', 'LineWidth', 2);
line(x_ends, (trial_idx+0.5)*[1 1], 'Color', 'w', 'LineWidth', 2);
title('All trials');

set(gui.raster, 'ButtonDownFcn', @jump_to_trial);

% Show trace
%------------------------------------------------------------
gui.trace = subplot(3,2,2);
t = trial_frames(1):trial_frames(end);
plot(t, trace, '.-', 'LineWidth', 1, 'HitTest', 'off');
xlim(t([1 end]));
ylim(trace_scale);
% grid on;
title(sprintf('Trial %d', trial_idx));
xlabel('Frames relative to gate close');
ylabel('Signal [a.u.]');
hold on;
plot(trial_frames(2)*[1 1], trace_scale, 'k--', 'HitTest', 'off'); % Open-gate
plot(trial_frames(3)*[1 1], trace_scale, 'k--', 'HitTest', 'off'); % Close-gate
gui.trace_bar = plot(trial_frames(1)*[1 1], trace_scale, 'r', 'ButtonDownFcn', @start_drag); % Vertical bar -- can be dragged
gui.trace_dot = plot(trial_frames(1), trace(1), 'r.',... % Dot
         'MarkerSize', 24,...
         'HitTest', 'off');

% Show events, if present
if ds.is_eventdata_loaded
    eventdata = ds.get_events(cell_idx, trial_idx);
    num_events = size(eventdata,1);
    if (num_events > 0)
        peak_times = t(eventdata(:,2));
        X = kron(peak_times, [1 1 NaN]);

        % Draw event locations
        Y = repmat([trace_scale NaN], 1, num_events);
        plot(X, Y, 'm:', 'HitTest', 'off');

        % Draw amplitudes
        Y = zeros(3, num_events);
        Y(1,:) = trace(eventdata(:,2));
        Y(2,:) = Y(1,:) - eventdata(:,3)';
        Y(3,:) = NaN;
        plot(X, Y(:), 'm', 'LineWidth', 2, 'HitTest', 'off');
    end
else
    num_events = 0;
end
hold off;

set(h_fig, 'WindowButtonUpFcn', @end_drag);
set(gui.trace, 'ButtonDownFcn', @update_frame);
     
% Show behavior movie
%------------------------------------------------------------
gui.behavior = subplot(3,2,[4 6]);
gui.behavior_image = imagesc(Mb(:,:,1));
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
    gui.behavior_pos = plot(centroids(1,1), centroids(1,2), 'ro');
    
    if (num_events > 0)
        peak_frames = eventdata(:,2);
        plot(centroids(peak_frames,1), centroids(peak_frames,2), 'm.', 'MarkerSize', 16);
    end
end

    function start_drag(~, ~)
        set(h_fig, 'WindowButtonMotionFcn', @update_frame);
    end

    function end_drag(~, ~)
        set(h_fig, 'WindowButtonMotionFcn', '');
    end

    function update_frame(~, e)
        sel_frame = round(e.IntersectionPoint(1));
        
        sel_frame = sel_frame + alignment_frame;
        sel_frame = max(sel_frame, 1);
        sel_frame = min(sel_frame, num_frames_in_trial);
        
        % Update visuals
        set(gui.trace_bar, 'XData', (sel_frame-alignment_frame)*[1 1]);
        set(gui.trace_dot, 'XData', (sel_frame-alignment_frame), 'YData', trace(sel_frame));
        set(gui.behavior_image, 'CData', Mb(:,:,sel_frame));
        
        subplot(gui.behavior);
        title(sprintf('Frame %d of %d', sel_frame, num_frames_in_trial));
        
        if (ds.is_tracking_loaded)
            centroid = ds.trials(trial_idx).centroids(sel_frame, :);
            set(gui.behavior_pos, 'XData', centroid(1), 'YData', centroid(2));
        end
    end % update_frame

    function jump_to_trial(~, e)
        trial_idx2 = round(e.IntersectionPoint(2));
        trial_idx2 = max(1, trial_idx2);
        trial_idx2 = min(trial_idx2, ds.num_trials);
        
        if (trial_idx2 ~= trial_idx)
            go_to_trial(trial_idx2);
        end
    end % jump_to_trial

% Navigation controls
%------------------------------------------------------------
back_btn = uicontrol('Style', 'pushbutton',...
    'String', '<<',...
    'Units', 'normalized',...
    'Position', [0 0.1 0.05 0.8],...
    'Callback', @back_to_raster);

prev_trial_btn = uicontrol('Style', 'pushbutton',...
    'String', '^^',...
    'Units', 'normalized',...
    'Position', [0 0.9 0.05 0.1],...
    'Callback', {@trial_btn_callback, trial_idx-1});
if trial_idx == 1
    prev_trial_btn.Enable = 'off';
end

next_trial_btn = uicontrol('Style', 'pushbutton',...
    'String', 'vv',...
    'Units', 'normalized',...
    'Position', [0 0.0 0.05 0.1],...
    'Callback', {@trial_btn_callback, trial_idx+1});
if trial_idx == ds.num_trials
    next_trial_btn.Enable = 'off';
end

    function back_to_raster(~, ~)
        % Clear callbacks associated with figure
        set(h_fig, 'WindowButtonUpFcn', '');
        view_raster_touch(ds, cell_idx, 'fig', h_fig, 'type', subraster_type);
    end

    function trial_btn_callback(~, ~, trial_idx)
        go_to_trial(trial_idx);
    end

    function go_to_trial(trial_idx)
        set(h_fig, 'WindowButtonUpFcn', '');
        view_detailed_trial_touch(ds, cell_idx, trial_idx, 'fig', h_fig, 'type', subraster_type);
    end

end % view_detailed_trial_touch

function [trial_frames, alignment_frame] = get_frames_relative_to_close(ds, trial_idx)
    trial_frames = double(ds.trial_indices(trial_idx, :)); % [start open-gate close-gate end]
    trial_frames = trial_frames - trial_frames(1) + 1;
    alignment_frame = trial_frames(3); % Align to close of gate

    trial_frames = trial_frames - alignment_frame;
end % get_frames_relative_to_close

