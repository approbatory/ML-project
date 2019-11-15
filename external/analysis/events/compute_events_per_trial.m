function events_per_trial = compute_events_per_trial(eventdata, trial_indices, full_num_frames)
% Reorganize event data by trials -- as required by DaySummary event
% storage. Namely,
%
% Given:
%   eventdata: {num_cells x 1} cell where,
%       eventdata{k}: [num_events x 3] is the event data for cell k
%
% Return:
%   events_per_trial: {num_trials x 1} cell where,
%       events_per_trial{t}: {num_cells x 1} contains the event data for
%       all cells in trial t. Frame numbers are referenced to each trial.
%
if nargin < 3
    full_num_frames = trial_indices(end, end);
end

num_trials = size(trial_indices, 1);

% A lookup table from frames to trial index
frames2trial = zeros(full_num_frames, 1);
for k = 1:num_trials
    frames2trial(trial_indices(k,1):trial_indices(k,end)) = k;
end

num_cells = length(eventdata);
events = cell(num_cells, num_trials);

for c = 1:num_cells
    eventdata_c = eventdata{c}; % Format: [trough-frame peak-frame amp]
    if ~isempty(eventdata_c)
        num_events = size(eventdata_c,1);
        events2trial = frames2trial(eventdata_c(:,2));

        % A trial "chunk" is the set of events in eventdata_c that belong
        % to the same trial. Note that we are assuming that the events in
        % eventdata_c are ordered in time.
        trial_chunk_start_inds = [1; find(diff(events2trial))+1];
        num_trial_chunks = length(trial_chunk_start_inds);
        for k = 1:num_trial_chunks
            % Pull out a single trial chunk from eventdata
            start_idx = trial_chunk_start_inds(k);
            if k ~= num_trial_chunks
                end_idx = trial_chunk_start_inds(k+1)-1;
            else
                end_idx = num_events;
            end
            trial_chunk = eventdata_c(start_idx:end_idx, :);

            % Note: A trial_idx of 0 indicates that the frames associated
            % with these events do not occur in 'trial_indices' -- which
            % can occur if 'trial_indices' has probe trials removed.
            trial_idx = events2trial(start_idx);
            if trial_idx ~= 0                
                % Reindex the frame indices on a per-trial basis
                trial_init_frame = trial_indices(trial_idx, 1);
                trial_chunk(:,1:2) = trial_chunk(:,1:2) - trial_init_frame + 1;

                events{c, trial_idx} = trial_chunk;
            end
        end
    end
end

events_per_trial = cell(num_trials, 1);
for k = 1:num_trials
    events_per_trial{k} = events(:,k);
end

end % compute_eventdata_per_trial