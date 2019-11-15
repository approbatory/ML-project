function plot_trace(obj, cell_idx)
    % Plot the trace of a single cell, color-coded by trial
    trace_min = Inf;
    trace_max = -Inf;

    colors = 'kbr';
    for k = 1:obj.num_trials
        [trace, inds] = obj.get_trace(cell_idx, k);

        trial_trace_min = min(trace);
        trial_trace_max = max(trace);
        if (trial_trace_min < trace_min)
            trace_min = trial_trace_min;
        end
        if (trial_trace_max > trace_max)
            trace_max = trial_trace_max;
        end

        plot(inds, trace,...
             colors(mod(k,length(colors))+1));
        hold on;
    end
    hold off;
    xlim([1 inds(end)]); % Using the last 'inds' from loop!
    trace_delta = trace_max - trace_min;
    ylim([trace_min trace_max] + 0.1*trace_delta*[-1 1]);
    xlabel('Frame index');
    ylabel('Signal [a.u.]');
end