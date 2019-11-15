function [global_min, global_max] = get_trace_bounds(ds, cell_idx)
    % Scan through all trials for a single cell trace from DaySummary,
    % so that we can apply a common scaling to the trace
    global_min = Inf; % Global min
    global_max = -Inf; % Global max
    
    for trial_idx = 1:ds.num_trials
        trace = ds.get_trace(cell_idx, trial_idx);
        m = min(trace);
        if m < global_min
            global_min = m;
        end
        M = max(trace);
        if M > global_max
            global_max = M;
        end
    end
end