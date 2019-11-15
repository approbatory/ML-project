function [resp, state] = view_cell_interactively(ds, cell_idx, movie, fps, state)
% Visually inspect the active portions of a trace side-by-side with
%   the provided miniscope movie.

filter = ds.cells(cell_idx).im;

[trace_orig, frames_to_movie] = ds.get_trace(cell_idx);
trace_fixed = fix_baseline(trace_orig);
time = 1/fps*((1:length(trace_orig))-1);

% Some parameters
filter_threshold = 0.3; % For generating the filter outline
num_neighbors_to_show = min(25, ds.num_cells-1);

% Display parameters
active_frame_padding = round(5*fps); % Used by 'parse_active_frames'
time_window = 100/fps; % Width of running window

% Set up GUI
%------------------------------------------------------------
set(state.fig_handle, 'WindowButtonUpFcn', @end_drag);
global_trace = subplot(3,3,[1 2 3]);
running_trace = subplot(3,3,[6 9]);
movie_subplot = subplot(3,3,[4 5 7 8]);

movie_clim = state.movie_clim;
h = imagesc(rescale_filter_to_clim(filter, movie_clim), movie_clim);
set(h, 'ButtonDownFcn', @add_point_of_interest);
colormap gray;
axis image;
xlabel('x [px]');
ylabel('y [px]');
hold on;

% Show boundary of current filter
boundaries = compute_ic_boundary(filter, filter_threshold);
for i = 1:length(boundaries)
    boundary = boundaries{i};
    plot(boundary(:,1), boundary(:,2), 'c', 'LineWidth', 2,...
         'HitTest', 'off');
end

% Plot boundaries of other cells, and retrieve their handles so that
%   we can toggle the boundaries on and off
%------------------------------------------------------------
% Old method plots ALL other cells
% other_cells = setdiff(1:ds.num_cells, cell_idx);
% num_other_cells = ds.num_cells - 1;

% New method only plots N nearest sources -- improves performance when
% there are 1000+ sources to classify
num_other_cells = min(300, ds.num_cells-1);
other_cells = ds.get_nearest_sources(cell_idx, num_other_cells);
other_cell_handles = zeros(num_other_cells, 1);

for n = 1:num_other_cells
    oc_idx = other_cells(n);
    boundary = ds.cells(oc_idx).boundary;
    if isempty(ds.cells(oc_idx).label)
        color = 'w';
    else
        if ds.is_cell(oc_idx)
            color = 'g';
        else
            color = 'r';
        end
    end
    other_cell_handles(n) = plot(boundary(:,1), boundary(:,2), color,...
                                 'HitTest', 'off');
end
show_map(state.show_map);

% Plot boundaries and indices of the nearest neighbors of the
% current cell
%------------------------------------------------------------
corr_colors = redblue(201);

neighbor_indices = ds.get_nearest_sources(cell_idx, num_neighbors_to_show);
neighbor_handles = zeros(num_neighbors_to_show, 2); % [Boundary Text]
for n = 1:num_neighbors_to_show
    neighbor_idx = neighbor_indices(n);
    boundary = ds.cells(neighbor_idx).boundary;
    com = ds.cells(neighbor_idx).com;
    
    % Compute the color from trace correlations
    trace_corr = ds.trace_corrs(cell_idx, neighbor_idx);
    trace_corr = round(trace_corr*100) + 101;
    color = corr_colors(trace_corr,:);
    
    neighbor_handles(n,1) = plot(boundary(:,1), boundary(:,2),...
                                 'HitTest', 'off',...
                                 'Clipping', 'on',...
                                 'Color', color,...
                                 'LineWidth', 2,...
                                 'HitTest', 'off');
    neighbor_handles(n,2) = text(com(1), com(2), num2str(neighbor_idx),...
                                 'HitTest', 'off',...
                                 'Clipping', 'on',...
                                 'HorizontalAlignment', 'center',...
                                 'Color', 'w',...
                                 'FontWeight', 'bold');
end
show_neighbors(state.show_neighbors);

% Plot points of interest
%------------------------------------------------------------
num_points = size(state.points_of_interest, 1);
for i = 1:num_points
    pt = state.points_of_interest(i,:);
    plot(pt(1), pt(2), 'y*');
end

% Indicate the center of mass of the current filter
COM = ds.cells(cell_idx).com;
plot(COM(1), COM(2), 'b.');
hold off;

% Start off zoomed
[height, width, ~] = size(movie);
zoom_half_width = min([width, height])/10;
xlim(COM(1)+zoom_half_width*[-1 1]);
ylim(COM(2)+zoom_half_width*[-1 1]);

% Compute the active portions of the trace and display
%------------------------------------------------------------
if state.baseline_removed
    trace = trace_fixed;
else
    trace = trace_orig;
end

thresh = compute_threshold(trace, state.threshold_scale);
[active_periods, num_active_periods] = ...
    parse_active_frames(trace > thresh, active_frame_padding);

setup_traces();

% Interaction loop:
%   Display the user-specified active period
%------------------------------------------------------------
prompt = 'Cell viewer >> ';
resp = lower(strtrim(input(prompt, 's')));
val = str2double(resp);

% State of interaction loop (will not carry over to other cells)
temp_state.last_val = [];
temp_state.zoomed = true;

while (1)
    if (~isnan(val)) % Is a number
        if ((1 <= val) && (val <= num_active_periods))
            display_active_period(val);
            temp_state.last_val = val;
        else
            fprintf('  Sorry, %d is not a valid period index for this trace\n', val);
        end
    else % Not a number
        switch (resp)
            case {'', 'q'} % "quit"
                % Unset the handler on fig handle, which persists!
                set(state.fig_handle, 'WindowScrollWheelFcn', '');
                break;
                
            case 't' % "threshold"
                fprintf('  Please select a new threshold on the global trace\n');
                while (1)
                    [~, thresh] = ginput(1);
                    if (gca == global_trace)
                        break;
                    else
                        fprintf('  Error! New threshold must be defined on the GLOBAL trace\n');
                    end
                end
                thresh_scale = (thresh - min(trace)) / (max(trace) - min(trace));
                fprintf('  New threshold value of %.3f selected (%.1f%% of max)!\n',...
                    thresh, thresh_scale*100);

                % Recompute active periods and redraw
                [active_periods, num_active_periods] =...
                    parse_active_frames(trace > thresh, active_frame_padding);
                setup_traces();
                
                % Persist the threshold
                if (0 < thresh_scale) && (thresh_scale < 1)
                    state.threshold_scale = thresh_scale;
                end
                
                temp_state.last_val = []; % Last displayed segment no longer valid
                
            case 'b' % Fix "baseline"
                state.baseline_removed = ~state.baseline_removed; % Toggle
                if (state.baseline_removed)
                    trace = trace_fixed;
                    fprintf('  Showing trace with baseline correction\n');
                else
                    trace = trace_orig;
                    fprintf('  Showing original trace without baseline correction\n');
                end
                
                thresh = compute_threshold(trace, state.threshold_scale);
                [active_periods, num_active_periods] = ...
                    parse_active_frames(trace > thresh, active_frame_padding);
                setup_traces();
                
                temp_state.last_val = []; % Last displayed segment no longer valid
                
            case 'r' % "replay"
                if ~isempty(temp_state.last_val)
                    display_active_period(temp_state.last_val);
                end
                            
            case 'z' % "zoom"
                subplot(movie_subplot); % Focus on the movie subplot
                if (temp_state.zoomed) % Return to original view
                    xlim([1 width]);
                    ylim([1 height]);
                    temp_state.zoomed = false;
                else
                    xlim(COM(1)+zoom_half_width*[-1 1]);
                    ylim(COM(2)+zoom_half_width*[-1 1]);
                    temp_state.zoomed = true;
                end
                    
            case {'h', 'l'} % "higher/lower contrast"
                subplot(movie_subplot); % Focus on the movie subplot
                c_range = diff(movie_clim);
                if (strcmp(resp, 'h'))
                    movie_clim = movie_clim + c_range*[0.1 -0.1];
                    fprintf('  Increased contrast (new CLim=[%.3f %.3f])\n',...
                        movie_clim(1), movie_clim(2));
                else
                    movie_clim = movie_clim + c_range*[-0.1 0.1];
                    fprintf('  Decreased contrast (new CLim=[%.3f %.3f])\n',...
                        movie_clim(1), movie_clim(2));
                end
                set(gca, 'CLim', movie_clim);
                state.movie_clim = movie_clim;
            
            case 'f' % Show "filter"
                set(h, 'CData', rescale_filter_to_clim(filter, state.movie_clim));
                
            case 'm' % Show "map" (i.e. all other cells)
                state.show_map = ~state.show_map;
                show_map(state.show_map);

            case 'x' % Show cross correlation of neighbors
                state.show_neighbors = ~state.show_neighbors;
                show_neighbors(state.show_neighbors);
                
            otherwise
                fprintf('  Sorry, could not parse "%s"\n', resp);
        end
    end

    resp = lower(strtrim(input(prompt, 's')));
    val = str2double(resp);
end

    % Display subroutines
    %------------------------------------------------------------
    function setup_traces()
        global t_g t_r dot;
        
        x_range = [time(1) time(end)];
        y_range = [min(trace(:)) max(trace(:))];
        y_delta = y_range(2) - y_range(1);
        y_range = y_range + 0.1*y_delta*[-1 1];
        
        % Prepare global trace
        subplot(global_trace)
        plot(time, trace, 'b', 'HitTest', 'off');
        hold on;
        plot(x_range, thresh*[1 1], 'r--', 'HitTest', 'off');
        for period_idx = 1:num_active_periods
            active_period = active_periods(period_idx, :);
            active_frames = active_period(1):active_period(2);
            plot(time(active_frames), trace(active_frames), 'r',...
                 'HitTest', 'off');
            text(double(time(active_frames(1))),... % 'text' fails on single
                 double(y_range(2)),...
                 num2str(period_idx),...
                 'Color', 'r',...
                 'VerticalAlignment', 'top',...
                 'HitTest', 'off');
        end
        xlim(x_range);
        ylim(y_range);
        t_g = plot(time(1)*[1 1], y_range, 'k', 'ButtonDownFcn', @start_drag); % Time indicator
        xlabel('Time [s]');
        ylabel('Signal [a.u.]');
        title(sprintf('Source %d of %d', cell_idx, ds.num_cells));
        hold off;
        
        % Prepare running trace
        subplot(running_trace);
        plot(time, trace, 'b', 'HitTest', 'off');
        hold on;
        for period_idx = 1:num_active_periods
            active_period = active_periods(period_idx, :);
            active_frames = active_period(1):active_period(2);
            plot(time(active_frames), trace(active_frames), 'r',...
                 'HitTest', 'off');
        end
        xlim([0 time_window]);
        ylim(y_range);
        t_r = plot(time(1)*[1 1], y_range, 'k', 'HitTest', 'off'); % Time indicator
        dot = plot(time(1), trace(1), 'or',...
                    'MarkerFaceColor', 'r',...
                    'MarkerSize', 12,...
                    'HitTest', 'off'); % Dot
        xlabel('Time [s]');
        ylabel('Signal [a.u.]');
        hold off;
        
        render_frame(1);
        
        % Event handlers for trace windows
        set(global_trace, 'ButtonDownFcn', @go_to_selected_frame);
        set(running_trace, 'ButtonDownFcn', @go_to_selected_frame);
        set(state.fig_handle, 'WindowScrollWheelFcn', @scroll_frame);
    end % setup_traces
    
    function display_active_period(selected_indices)
        global break_active_period;
        
        % The following flag allows for premature termination of active
        % period display. Flag is set by 'go_to_selected_frame'
        break_active_period = false;
        for selected_idx = selected_indices
            frames = active_periods(selected_idx,1):...
                     active_periods(selected_idx,2);
            for k = frames
                if break_active_period
                    break; %#ok<UNRCH>
                end
                render_frame(k);
            end
        end
    end % display_active_period

    function show_map(show)
        vis_val = 'off';
        if (show)
            vis_val = 'on';
        end
        
        for m = 1:num_other_cells
            set(other_cell_handles(m), 'Visible', vis_val);
        end
    end % show_map

    function show_neighbors(show)
        vis_val = 'off';
        if (show)
            vis_val = 'on';
        end
        
        for m = 1:num_neighbors_to_show
            set(neighbor_handles(m,1), 'Visible', vis_val);
            set(neighbor_handles(m,2), 'Visible', vis_val);
        end
    end

    function add_point_of_interest(~, e)
        coord = round(e.IntersectionPoint([1 2]));
        
        % Add new point of interest to plot
        subplot(movie_subplot);
        hold on;
        plot(coord(1), coord(2), 'y*');
        hold off;
        
        frame_idx = get(h, 'UserData');
        state.points_of_interest = [state.points_of_interest; coord frame_idx];
    end

    function go_to_selected_frame(~, e)
        global break_active_period;
        break_active_period = true;
        t = e.IntersectionPoint(1);
        k = round(fps*t+1);
        render_frame(k);
    end

    function start_drag(~,~)
        set(state.fig_handle, 'WindowButtonMotionFcn', @drag_frame);
    end

    function end_drag(~,~)
        set(state.fig_handle, 'WindowButtonMotionFcn', '');
    end

    function drag_frame(~,~)
        cp = get(global_trace, 'CurrentPoint');
        t = cp(1);
        k = round(fps*t+1);
        render_frame(k);
    end

    function render_frame(k)
        global t_g t_r dot;
        k = max(1,k); k = min(length(frames_to_movie),k); % Clamp
        A = movie(:,:,frames_to_movie(k));
        set(h, 'CData', A);
        set(h, 'UserData', k); % So that we can refer to the displayed frame elsewhere
        set(t_g, 'XData', time(k)*[1 1]);
        set(t_r, 'XData', time(k)*[1 1]);
        set(dot, 'XData', time(k), 'YData', trace(k));
        set(running_trace, 'XLim', time(k) + time_window/2*[-1 1]);
        drawnow;
    end

    function scroll_frame(~, e)
        current_frame = get(h, 'UserData');
        if (e.VerticalScrollCount < 0) % Scroll up
            k = current_frame - 1;
        else
            k = current_frame + 1;
        end
        render_frame(k);
    end
end % main function

function [active_frames, num_active] = parse_active_frames(binary_trace, half_width)
% Segment the active portions of a binary trace into intervals

    if (half_width > 0)
        trace = single(binary_trace);
        trace = conv(trace, ones(1, 2*half_width + 1), 'same');
        trace = logical(trace);
    else
        trace = binary_trace;
    end
    trace_comp = ~trace; % Complement of the trace

    active_frames = [];

    % Loop to find all activity transitions in the trace
    curr = 1;
    while (1)
        next = find(trace(curr:end), 1, 'first');
        if (isempty(next))
            break;
        end
        active_frames = [active_frames curr+(next-1)]; %#ok<*AGROW>
        curr = curr + next;

        next = find(trace_comp(curr:end), 1, 'first');
        if (isempty(next))
            break;
        end
        active_frames = [active_frames curr+(next-2)];
        curr = curr + next;
    end

    if (mod(length(active_frames),2) == 1) % Ended with active frame
        active_frames = [active_frames length(binary_trace)];
    end

    active_frames = reshape(active_frames, 2, length(active_frames)/2)';
    num_active = size(active_frames, 1);
end

function thresh = compute_threshold(tr, scale)
    trace_max = max(tr);
    trace_min = min(tr);
    thresh = scale * (trace_max - trace_min) + trace_min;
end

function filter_out = rescale_filter_to_clim(filter, clim)
% Numerically rescale the filter to match the provided clim

    f_max = max(filter(:));
    f_min = min(filter(:));

    filter_norm = (filter-f_min)/f_max; % Matched to range [0, 1]

    clim_delta = clim(2)-clim(1);
    clim_usage = 0.8;
    filter_out = clim(1) + (1-clim_usage)/2*clim_delta +...
                 clim_usage*clim_delta*filter_norm;

end
