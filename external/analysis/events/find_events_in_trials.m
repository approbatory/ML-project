function events = find_events_in_trials(trace, trial_indices, threshold, baseline, amp_threshold)
% Small wrapper around 'find_events' so that event parameter determination
% (e.g. trough preceding the event peak, event amplitude) respects trial
% boundaries in the trace. This is necessary in the prefrontal dataset,
% since the traces are concatenations of multiple trials -- time is not
% actually continuous across the trial boundaries.
%

if (nargin < 5)
    amp_threshold = 0;
end

num_trials = size(trial_indices, 1);

events = cell(num_trials, 1);
for k = 1:num_trials
    % Find events within the trial
    trial_frames = trial_indices(k,[1 end]);
    trial_frames = trial_frames(1):trial_frames(2);
    es = find_events(trace(trial_frames), threshold, baseline);
    
    % Add frame offset for the k-th trial
    if ~isempty(es)
        es(:,1) = es(:,1) + trial_frames(1) - 1;
        es(:,2) = es(:,2) + trial_frames(1) - 1;
    end
    events{k} = es;
end

events = cell2mat(events);

% Filter for amplitude heights
if ~isempty(events)
    max_event_amplitude = max(events(:,3));
    filtered_events = events(:,3) > amp_threshold * max_event_amplitude;
    events = events(filtered_events,:);
end