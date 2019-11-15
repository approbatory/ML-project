function plot_superposed_trials(obj, cell_idx, varargin)
    % Optional arguments allow for filtering of trials, e.g.
    %   "plot_superposed_trials(cell_idx, 'start', 'east')"
    display_trial = ones(obj.num_trials, 1);
    if ~isempty(varargin)
        display_trial = obj.filter_trials(varargin{:});
    end

    trace_min = Inf;
    trace_max = -Inf;

    colors = 'kbr';
    for k = 1:obj.num_trials
        if display_trial(k)
            trial_trace = obj.get_trace(cell_idx, k);

            trial_trace_min = min(trial_trace);
            trial_trace_max = max(trial_trace);
            if (trial_trace_min < trace_min)
                trace_min = trial_trace_min;
            end
            if (trial_trace_max > trace_max)
                trace_max = trial_trace_max;
            end

            plot(linspace(0, 1, length(trial_trace)),...
                 trial_trace,...
                 colors(mod(k,length(colors))+1));
            hold on;
        end
    end
    hold off;
    grid on;
    xlim([0 1]);
    trace_delta = trace_max - trace_min;
    ylim([trace_min trace_max] + 0.1*trace_delta*[-1 1]);
    xlabel('Trial phase [a.u.]');
    ylabel('Trace [a.u.]');
end