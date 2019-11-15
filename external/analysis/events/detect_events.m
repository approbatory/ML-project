function events = detect_events(ds, cell_idx, varargin)
% For formatting info on event data, see 'find_events.m'

use_prompt = true;
use_filter = true;
M = [];
movie_clim = [];
fps = 10;
cutoff_freq = [];
ext_events = [];

for j = 1:length(varargin)
    vararg = varargin{j};
    if ischar(vararg)
        switch lower(vararg)
            case 'noprompt'
                use_prompt = false;
            case 'fps'
                fps = varargin{j+1};
            case 'cutoff'
                cutoff_freq = varargin{j+1};
            case 'nofilter'
                use_filter = false;
            case {'m', 'movie'}
                M = varargin{j+1};
            case {'clim', 'movie_clim'}
                movie_clim = varargin{j+1};
            case {'ext', 'extevents', 'events'}
                ext_events = varargin{j+1};
        end
    end
end

trace_orig = ds.get_trace(cell_idx);
if ~isempty(M) && isempty(movie_clim)
    movie_clim = compute_movie_scale(M);
end

% We'll for events in a smoothed version of the trace
% Default parameters comes from cerebellar processing, where we used
%   - 30 Hz sampling frequency
%   - 4 Hz cutoff frequency
if isempty(cutoff_freq)
    cutoff_freq = 4/30 * fps;
end
if use_filter
%     fprintf('Applying LPF (fc=%.1f Hz) to trace...\n', cutoff_freq);

    % Don't apply LPF across trial boundaries
    filt_traces = cell(1,ds.num_trials);
    for tidx = 1:ds.num_trials
        filt_traces{tidx} = filter_trace(ds.get_trace(cell_idx,tidx), cutoff_freq, fps);
    end
    trace = cell2mat(filt_traces);
else
    trace = trace_orig;
end

% Basic trace properties
%------------------------------------------------------------
num_frames = length(trace);
[baseline, sigma, stats] = estimate_baseline_sigma(trace);

% Application state
state.allow_manual_events = false;
state.x_anchor = 1;
state.x_range = min(500, num_frames);
state.show_orig = use_prompt;%true;
state.show_dots = false;
state.show_trials = use_prompt && (ds.num_trials > 1);
state.sel_event = 0;
state.last_requested_trial = 0;

init_info = struct('num_frames', num_frames,...
                   'baseline', baseline,...
                   'sigma', sigma,...
                   'threshold', baseline + 5*sigma,...
                   'amp_threshold', 0.1);
events = struct('info', init_info, 'auto', [], 'manual', []);

% FIXME: It'd be nice to not draw the figure when 'use_prompt' is disabled
hfig = figure;
if ~use_prompt
    set(hfig,'Visible','off');
end
gui = setup_gui(hfig, num_frames, compute_display_range(trace), stats, trace_orig);

if ~isempty(ext_events)
    if isstruct(ext_events)
        % If struct, assume to have originated from the "default" event
        % detection method implemented here and in 'find_events_in_trials'
        init_info = ext_events.info;
    else        
        % TODO: Handle different types of external event specification
        events.manual = ext_events;
        redraw_manual_events(gui);
    end
end

set_threshold([], [], gui);
update_gui_state(gui, state);
if state.show_trials
    set_trial(1, gui);
end

% Interaction loop:
%------------------------------------------------------------
prompt = 'Detector >> ';
while (use_prompt)
    resp = lower(strtrim(input(prompt, 's')));
    val = str2double(resp);

    if (~isnan(val)) % Is a number
        if (1 <= val) && (val <= ds.num_trials)
            set_trial(val, gui);
        end
    else % Not a number
        switch (resp)
            case 'q' % "quit"
                break;

            case 'z' % zoom in
                x_center = state.x_anchor + 1/2 * state.x_range;
                state.x_range = 0.5*state.x_range;
                state.x_anchor = x_center - 1/2 * state.x_range;
                redraw_local_window(gui, state);
                
            case 'o' % zoom out
                x_center = state.x_anchor + 1/2 * state.x_range;
                state.x_range = 2*state.x_range;
                state.x_anchor = x_center - 1/2 * state.x_range;
                redraw_local_window(gui, state);
                
            case 'r' % toggle "raw"/original trace
                state.show_orig = ~state.show_orig;
                update_gui_state(gui, state);
                
            case 'd' % show dots
                state.show_dots = ~state.show_dots;
                update_gui_state(gui, state);
                
            case 'b' % show trials
                state.show_trials = ~state.show_trials;
                update_gui_state(gui, state);
                
            case {'', 'n'} % next trial
                if state.last_requested_trial < ds.num_trials
                    set_trial(state.last_requested_trial+1, gui);
                end
                
            case 't' % reset threshold
                events.info = init_info;
                set_threshold([], [], gui);
                                
%             case 'm' % toggle manual input -- DISABLED for now
%                 state.allow_manual_events = ~state.allow_manual_events;
                
%             case 'x' % erase last event
%                 num_events = length(events.manual);
%                 if (num_events > 0)
%                     events.manual = events.manual(1:num_events-1);
%                     redraw_manual_events(gui);
%                 else
%                     fprintf('  No events to remove!\n');
%                 end

            otherwise
                fprintf('  Sorry, could not parse "%s"\n', resp);
        end
    end
end % Main interaction loop

% Finished interaction -- prepare output
%------------------------------------------------------------
close(hfig);

% Re-sort auto events by peak event time
if ~isempty(events.auto)
    events.auto = sortrows(events.auto, 2);
end

    % Supplementary functions
    %------------------------------------------------------------
    function get_next_page(gui)
        current_end = state.x_anchor + state.x_range;
        if (current_end >= gui.num_frames)
            new_anchor = gui.num_frames - state.x_range + 1;
        else
            new_anchor = state.x_anchor + 0.1*state.x_range + 1;
        end

        state.x_anchor = new_anchor;
        redraw_local_window(gui, state);
    end % get_next_page

    function get_prev_page(gui)
        new_anchor = state.x_anchor - (0.1*state.x_range + 1);
        state.x_anchor = max(1, new_anchor);
        redraw_local_window(gui, state);
    end % get_prev_page

    function gui = setup_gui(hf, num_frames, trace_display_range, stats, trace_orig)
        % Display parameters kept around for convenience
        gui.hfig = hf;
        gui.num_frames = num_frames;
        gui.trace_display_range = trace_display_range;
        
        % Setup the GLOBAL trace plot
        %------------------------------------------------------------
        gui.global = subplot(2,4,1:3);
        gui.global_rect = rectangle('Position',[-1 trace_display_range(1) 1 diff(trace_display_range)],...
                  'EdgeColor', 'none',...
                  'FaceColor', 'c', 'HitTest', 'off');
        hold on;
        plot(trace, 'k', 'HitTest', 'off');
        gui.global_thresh = plot([1 num_frames], -Inf*[1 1], 'm--', 'HitTest', 'off');
        gui.global_auto = plot(-Inf, -1, 'm.', 'HitTest', 'off');
        gui.global_manual = plot(-Inf, -1, 'r.', 'HitTest', 'off');
        hold off;
        box on;
        xlim([1 num_frames]);
        ylim(trace_display_range);
        xlabel('Frame');
        ylabel('Fluorescence (a.u.)');
        
        % Setup the HISTOGRAM plot
        %------------------------------------------------------------
        gui.histogram = subplot(4,4,4);
        semilogy(stats.hist_centers, stats.hist_counts, 'k.', 'HitTest', 'off');
        xlim(trace_display_range);
        hold on;
        % First power of 10 that exceeds the maximum count
        count_range = [1 10^ceil(log10(max(stats.hist_counts)))];
        ylim(count_range);
        for k = 1:size(stats.percentiles,1)
            f = stats.percentiles(k,2);
            plot(f*[1 1], count_range, 'Color', 0.5*[1 1 1], 'HitTest', 'off');
        end
        gui.histogram_thresh = plot(-Inf*[1 1], count_range, 'm', 'HitTest', 'off');
        hold off;
        xlabel('Fluorescence');
        ylabel('Trace histogram');

        % Setup the event amplitude CDF
        %------------------------------------------------------------
        gui.event_amp_cdf = subplot(4,4,8);
        gui.cdf = plot(-1, -1, 'm.-', 'HitTest', 'off');
        hold on;
        gui.cdf_sel_event = plot(-1, -1, 'mo', 'HitTest', 'off');
        gui.cdf_amp_threshold = plot(-1*[1 1], [0 1], 'k', 'HitTest', 'off');
        hold off;
        xlim([0 1]);
        ylim([0 1]);
        xticks(0:0.1:1);
        yticks(0:0.1:1);
        grid on;
        xlabel('Norm event amplitude');
        ylabel('CDF');
        
        % Setup the LOCAL trace plot
        %------------------------------------------------------------
        if ~isempty(M) % Show movie
            gui.movie = subplot(2,4,8);
            gui.movie_frame = imagesc(ds.cells(cell_idx).im,...
                movie_clim);
            axis image;
            colormap gray;
            hold on;
            boundary = ds.cells(cell_idx).boundary;
            plot(boundary(:,1), boundary(:,2), 'c', 'LineWidth', 2, 'HitTest', 'off');
            com = ds.cells(cell_idx).com;
            plot(com(1), com(2), 'b.');
            num_neighbors = min(10, ds.num_classified_cells-1);
            neighbor_inds = ds.get_nearest_sources(cell_idx, num_neighbors);
            for n = neighbor_inds
                boundary = ds.cells(n).boundary;
                plot(boundary(:,1), boundary(:,2), 'w');
                
                ncom = ds.cells(n).com;
                text(ncom(1), ncom(2), num2str(n),...
                     'Color', 'w',...
                     'HorizontalAlignment', 'center',...
                     'HitTest', 'off', 'Clipping', 'on');
            end
            hold off;
            
            [height, width, ~] = size(M);
            zoom_half_width = min([height, width])/15;
            xlim(com(1) + zoom_half_width*[-1 1]);
            ylim(com(2) + zoom_half_width*[-1 1]);
            
            gui.local = subplot(2,4,5:7);
        else
            gui.local = subplot(2,1,2);
        end
        gui.local_orig = plot(trace_orig, 'Color', 0.6*[1 1 1], 'HitTest', 'off');
        hold on;
        plot(trace, 'k', 'HitTest', 'off');
        gui.local_dots = plot(trace, 'k.', 'HitTest', 'off');
        trial_starts = ds.trial_indices(:,1);
        gui.local_trials = plot(trial_starts, trace(trial_starts), 'ko', 'HitTest', 'off');
        text_y = trace_display_range(1) + 0.05*diff(trace_display_range);
        gui.local_trials_text = cell(ds.num_trials,1);
        for k = 1:ds.num_trials
            trial_start_k = double(trial_starts(k));
            gui.local_trials_text{k} = text(trial_start_k, text_y,...
                 sprintf('Trial %d', k),...
                 'HorizontalAlignment', 'center',...
                 'HitTest', 'off', 'Clipping', 'on');
        end
        gui.local_cursor_dot = plot(-Inf,trace(1),'ro',...
            'MarkerFaceColor','r',...
            'MarkerSize',6,'HitTest','off');
        gui.local_cursor_bar = plot(-Inf*[1 1], trace_display_range, 'k--', 'HitTest', 'off');
        gui.local_thresh = plot([1 num_frames], -Inf*[1 1], 'm--', 'HitTest', 'off');
        gui.local_auto = plot(-1, -1, 'm:', 'HitTest', 'off');
        gui.local_sel_event = rectangle('Position',[-1 trace_display_range(1) 0 diff(trace_display_range)],...
                  'EdgeColor', 'none',...
                  'FaceColor', [1 0 1 0.2],... // Transparent magenta
                  'HitTest', 'off');
        gui.local_sel_event_peak = plot([-1 -1], trace_display_range, 'm', 'HitTest', 'off');
        gui.local_auto_amps = plot(-1, -1, 'm', 'LineWidth', 2, 'HitTest', 'off');
        gui.local_manual = plot(-1, -1, 'r');
        hold off;
        ylim(trace_display_range);
        xlabel('Frame');
        ylabel('Fluorescence');
        
        % Add GUI event listeners
        set(gui.global, 'ButtonDownFcn', {@global_plot_handler, gui});
        set(gui.histogram, 'ButtonDownFcn', {@histogram_handler, gui});
        set(gui.event_amp_cdf, 'ButtonDownFcn', {@cdf_handler, gui});
        set(gui.local, 'ButtonDownFcn', {@local_plot_handler, gui});
        set(hf, 'WindowButtonMotionFcn', {@track_cursor, gui});
        set(hf, 'WindowScrollWheelFcn', {@scroll_plot, gui});
        
        function track_cursor(~, e, gui)
            x = round(e.IntersectionPoint(1));
            if ((1<=x)&&(x<=gui.num_frames))
                if state.allow_manual_events
                    x = seek_localmax(trace, x);
                end
                if ~isempty(M)
                    set(gui.movie_frame, 'CData', M(:,:,x));
                end
                set(gui.local_cursor_bar,'XData',x*[1 1]);
                set(gui.local_cursor_dot,'XData',x,'YData',trace(x));
            end
        end % track_cursor
        
        function scroll_plot(~, e, gui)
            if (state.show_trials && (ds.num_trials > 1)) % Scroll by trials
                trial_idx = state.last_requested_trial;
                if (e.VerticalScrollCount < 0) % Scroll up
                    trial_idx = trial_idx - 1;
                else
                    trial_idx = trial_idx + 1;
                end
                
                trial_idx = max(1, trial_idx); % Clamp
                trial_idx = min(trial_idx, ds.num_trials);
                set_trial(trial_idx, gui);
                
            else % Default scrolling
                if (e.VerticalScrollCount < 0) % Scroll up
                    get_prev_page(gui);
                else
                    get_next_page(gui);
                end
            end
        end % scroll_plot
        
    end % setup_gui

    % Update the GUI
    %------------------------------------------------------------
    function update_gui_state(gui, state)
        redraw_local_window(gui, state);
        
        if state.show_orig
            set(gui.local_orig, 'Visible', 'on');
        else
            set(gui.local_orig, 'Visible', 'off');
        end
        
        if state.show_dots
            set(gui.local_dots, 'Visible', 'on');
        else
            set(gui.local_dots, 'Visible', 'off');
        end
        
        if state.show_trials
            trials_vis = 'on';
        else
            trials_vis = 'off';
        end
        set(gui.local_trials, 'Visible', trials_vis);
        for k = 1:ds.num_trials
            set(gui.local_trials_text{k}, 'Visible', trials_vis);
        end
    end

    function redraw_local_window(gui, state)
        if use_prompt
            figure(gui.hfig);
        end
        rect_pos = get(gui.global_rect, 'Position');
        rect_pos(1) = state.x_anchor;
        rect_pos(3) = state.x_range;
        set(gui.global_rect, 'Position', rect_pos);
               
        subplot(gui.local);
        xlim([state.x_anchor, state.x_anchor+state.x_range-1]);
    end
    
    function redraw_threshold(gui)
        % Handle the possibility that the eventdata is empty
        if ~isempty(events.auto)
            auto_peaks = events.auto(:,2);
            num_auto_events = length(auto_peaks);
            cdf_x = events.auto(:,3) / events.auto(end,3); % Assume events are sorted by amplitude
            cdf_y = (1:num_auto_events)/num_auto_events;
        else
            auto_peaks = [];
            num_auto_events = 0;
            cdf_x = [];
            cdf_y = [];
        end
        
        % GLOBAL subplot
        set(gui.global_thresh, 'YData', events.info.threshold*[1 1]);
        set(gui.global_auto, 'XData', auto_peaks, 'YData', trace(auto_peaks));
        update_event_tally(gui);
        
        % HISTOGRAM subplot
        set(gui.histogram_thresh, 'XData', events.info.threshold*[1 1]);
        
        % CDF subplot
        set(gui.cdf_amp_threshold, 'XData', events.info.amp_threshold*[1 1]);
        set(gui.cdf, 'XData', cdf_x, 'YData', cdf_y);
        
        % LOCAL subplot
        set(gui.local_thresh, 'YData', events.info.threshold*[1 1]);
        
        % Note: NaN's break connections between line segments
        X = kron(auto_peaks', [1 1 NaN]);
        Y = repmat([gui.trace_display_range NaN], 1, num_auto_events);
        set(gui.local_auto, 'XData', X, 'YData', Y);
        
        % Draw event amplitudes
        if (num_auto_events > 0)
            Y = zeros(3, num_auto_events);
            Y(1,:) = trace(auto_peaks);
            Y(2,:) = Y(1,:) - events.auto(:,3)'; % Peak minus amplitude
            Y(3,:) = NaN;
            Y = Y(:);
        else
            Y = [];
        end
        set(gui.local_auto_amps, 'XData', X, 'YData', Y);
    end % redraw_threshold

    function redraw_manual_events(gui)
        set(gui.global_manual, 'XData', events.manual, 'YData', trace(events.manual));
        
        X = kron(events.manual', [1 1 NaN]);
        Y = repmat([gui.trace_display_range NaN], 1, length(events.manual));
        set(gui.local_manual, 'XData', X, 'YData', Y);
        
        update_event_tally(gui);
    end % redraw_manual_events

    function update_event_tally(gui)
        num_auto = size(events.auto,1);
%         num_manual = length(events.manual);
        
        subplot(gui.global);
        title(sprintf('Cell %d: %d events (auto only)', cell_idx, num_auto));
    end % update_event_tally

    % Event handlers for mouse input
    %------------------------------------------------------------
    function global_plot_handler(~, e, gui)
        switch e.Button
            case 1 % Left click -- Move the local viewpoint
                x = round(e.IntersectionPoint(1));
                if ((1<=x) && (x<=gui.num_frames))
                    if state.show_trials
                        t = find(x >= ds.trial_indices(:,1), 1, 'last');
                        set_trial(t, gui);
                    else
                        state.x_anchor = x - state.x_range/2;
                        redraw_local_window(gui, state);
                    end
                else
                    fprintf('\n  Not a valid frame for this trace!\n');
                end
                
            case 3 % Right click -- Set threshold
                t = e.IntersectionPoint(2);
                set_threshold(t, [], gui);
        end
    end % global_plot_handler

    function histogram_handler(~, e, gui)
        switch e.Button
            case 1 % Left click
                
            case 3 % Right click -- Set threshold
                t = e.IntersectionPoint(1);
                set_threshold(t, [], gui);
        end
    end % histogram_handler

    function cdf_handler(~, e, gui)
        switch e.Button
            case 1 % Left click -- Select a particular event
                x = e.IntersectionPoint(1);
                
                % Find the event with the nearest amplitude
                event_amps = events.auto(:,3);
                delta_amp = abs(event_amps - max(event_amps)*x);
                [~, se] = min(delta_amp);
                
                select_event(se, gui);
                
                % Move the local window
                sel_frame = events.auto(se,2);
                state.x_anchor = sel_frame - 1/2 * state.x_range;             
                redraw_local_window(gui, state);
            case 3 % Right click -- Set the amplitude threshold
                x = e.IntersectionPoint(1);
                x = max(0, x); x = min(x, 1); % Clamp to [0 1]
                
                set_threshold([], x, gui);
        end
    end % cdf_handler

    function local_plot_handler(~, e, gui)
        switch e.Button
            case 1 % Left click
                x = round(e.IntersectionPoint(1));
                if state.allow_manual_events
                    add_manual_event(x, gui);
                else
                    % Find the nearest event
                    event_times = events.auto(:,2);
                    delta_times = abs(event_times - x);
                    [~, se] = min(delta_times);
                    
                    select_event(se, gui);
                end
            case 3 % Right click
                
        end
    end % local_plot_handler

    % Data processing
    %------------------------------------------------------------
    function add_manual_event(x, gui)
        if ((1<=x) && (x<=gui.num_frames))
            x = seek_localmax(trace, x);
            % Don't make duplicate events
            auto_peaks = events.auto(:,2);
            if ~ismember(x, auto_peaks) && ~ismember(x, events.manual)
                events.manual = [events.manual; x];
                redraw_manual_events(gui);
            end
        else
            fprintf('\n  Not a valid event for this trace!\n');
        end
    end % add_manual_event

    function set_threshold(threshold, amp_threshold, gui)
        if ~isempty(threshold)
            events.info.threshold = threshold;
        end
        if ~isempty(amp_threshold)
            events.info.amp_threshold = amp_threshold;
        end

        i = events.info;
        events.auto = find_events_in_trials(trace, ds.trial_indices,...
            i.threshold, i.baseline, i.amp_threshold);
        if ~isempty(events.auto)
            events.auto = sortrows(events.auto, 3); % Sort events by amplitude
        end
        
        select_event(0, gui);
        redraw_threshold(gui);
    end % set_threshold

    function set_trial(trial_idx, gui)
        trial_start = ds.trial_indices(trial_idx,1);
        trial_end = ds.trial_indices(trial_idx,4);
        trial_range = trial_end - trial_start + 1;

        state.last_requested_trial = trial_idx;
        state.x_anchor = trial_start - 0.25 * trial_range;
        state.x_range = 1.5 * trial_range;
        redraw_local_window(gui, state);
    end % set_trial

    function select_event(event_idx, gui)
        % Event index refers to the row of 'events.auto'. An index of "0"
        % corresponds to not selecting any event.
        num_events = size(events.auto, 1);
        if (0 <= event_idx) && (event_idx <= num_events)
            state.sel_event = event_idx;
            if (event_idx > 0)
                event_amps = events.auto(:,3);
                sel_event_onset_frame = events.auto(event_idx,1);
                sel_event_peak_frame = events.auto(event_idx,2);
                sel_event_duration = sel_event_peak_frame - sel_event_onset_frame;
                sel_event_normamp = event_amps(event_idx)/max(event_amps);
                sel_event_cdf = event_idx / num_events;
            else
                sel_event_onset_frame = -1;
                sel_event_peak_frame = -Inf;
                sel_event_duration = 0;
                sel_event_normamp = -Inf;
                sel_event_cdf = -1;
            end
            set(gui.cdf_sel_event, 'XData', sel_event_normamp, 'YData', sel_event_cdf);
            rect_pos = get(gui.global_rect, 'Position');
            rect_pos(1) = sel_event_onset_frame;
            rect_pos(3) = sel_event_duration;
            set(gui.local_sel_event, 'Position', rect_pos);
            set(gui.local_sel_event_peak', 'Xdata', sel_event_peak_frame*[1 1]);
        end
    end % select_event

end % detect_events_interactively

function display_range = compute_display_range(trace)
    display_range = [min(trace) max(trace)];
    display_range = display_range + 0.1*diff(display_range)*[-1 1];
end % compute_display_range
