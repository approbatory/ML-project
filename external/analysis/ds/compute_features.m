function fv = compute_features(ds, cell_idx, trial_indices)
% Compute features from specified trials

num_trials = length(trial_indices);

fv = struct('times', [ds.trials(trial_indices).time]',...
            'mean_fluorescence', zeros(num_trials,1),...
            'max_fluorescence', zeros(num_trials,1),...
            'event_sum', zeros(num_trials,1),...
            'num_events', zeros(num_trials,1));
       
for k = 1:num_trials
    trial_idx = trial_indices(k);
    
    tr = ds.get_trace(cell_idx, trial_idx);
    
    % Basic features computed directly from fluorescence trace
    fv.mean_fluorescence(k) = mean(tr);
    fv.max_fluorescence(k) = max(tr);
    
    % Features based on detected events
    eventdata = ds.get_events(cell_idx, trial_idx);
    if ~isempty(eventdata)
        fv.event_sum(k) = sum(eventdata(:,3));
        fv.num_events(k) = size(eventdata, 1);
    end
end