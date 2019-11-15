function [raster, info] = plot_cell_raster(obj, cell_idx, varargin)
% Plots a raster of cell activity, where trials are aligned to the closing
% gate frame by default. Other alignment points are also possible via the 
% 'align' vararg. Additional arguments allow for filtering of trials, e.g.
%
%   plot_cell_raster(cell_idx, 'start', 'east')
%
% See 'DaySummary.filter_trials' for full details.
%
% Optional argument 'draw_correct' will place a box at the end
% of each trial indicating correct (green) or incorrect (red)
%
    draw_events = false;
    draw_correct = false;
    
    align_idx = 3; % By default, align to closing of gate   
    kept_trials = true(1, obj.num_trials);
    
    if ~isempty(varargin)
        for k = 1:length(varargin)
            vararg = varargin{k};
            if ischar(vararg)
                switch lower(vararg)
                    case {'align', 'align_to'}
                        align_idx = varargin{k+1};
                    case {'draw_events', 'event', 'events'}
                        draw_events = true;
                    case 'draw_correct'
                        draw_correct = true;
                end
            end
        end
        
        kept_trials = obj.filter_trials(varargin{:});
    end

    trial_inds = find(kept_trials);
    num_trials = length(trial_inds);
    if (num_trials == 0)
        % Display message in place of showing raster
        xlim([0 1]);
        ylim([0 1]);
        set(gca, 'XTick', []);
        set(gca, 'YTick', []);
        text(0.5, 0.5, {'No trials','selected'},...
             'HorizontalAlignment', 'center');
        return;
    end
    
    alignment_frames = obj.trial_indices(trial_inds, align_idx);
    [raster, info] = obj.get_aligned_trace(cell_idx, trial_inds, alignment_frames);
    
    switch align_idx
        case 1
            align_str = 'Frames relative to trial start';
        case 2
            align_str = 'Frames relative to gate open';
        case 3
            align_str = 'Frames relative to gate close';
        case 4
            align_str = 'Frames relative to trial end';
    end
    
    imagesc(info.aligned_time, 1:num_trials, raster, 'HitTest', 'off');
    set(gca, 'CLim', obj.trace_range(cell_idx,:));
    hold on;
    colormap parula; freezeColors;
    xlim(info.aligned_time([1 end]));
    xlabel(align_str);
    ylabel('Trial index');
    set(gca, 'TickLength', [0 0]);

    if (draw_events && obj.is_eventdata_loaded)
        for k = 1:num_trials
            trial_idx = trial_inds(k);
            eventdata = obj.get_events(cell_idx, trial_idx, 'align_to', align_idx);
            if ~isempty(eventdata)
                event_times = eventdata(:,2); % Peak times
                num_events = length(event_times);
                plot(event_times, k*ones(1,num_events), 'o',...
                     'MarkerSize', 4,...
                     'MarkerEdgeColor', 'k',...
                     'MarkerFaceColor', 'm',...
                     'HitTest', 'off');
            end
        end
    end
    
    if (draw_correct)
        corr_width = 0.025*size(raster,2);
        for k = 1:num_trials
            trial_idx = trial_inds(k);
            if obj.trials(trial_idx).correct
                corr_color = 'g';
            else
                corr_color = 'r';
            end
            rectangle('Position', [info.aligned_time(end) k-0.5 corr_width 1],...
                      'FaceColor', corr_color);
        end
        xlim([info.aligned_time(1) info.aligned_time(end)+corr_width]);
    end
    hold off;
end