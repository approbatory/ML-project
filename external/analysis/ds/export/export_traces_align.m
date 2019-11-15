function [X, xs, ys, alignment_axis] = export_traces_align(md, trial_map, align_idx)
% X = EXPORT_TRACES_ALIGN(md, trial_map, align_idx)
%
% Exports traces of all trials specified by trial_map into a 3D matrix 'X'
%   with dimensions [neurons x time x trials]. Trials are aligned by one of
%   four intra-trial events:
%     (1) Beginning of trial
%     (2) Opening of gate
%     (3) Closing of gate
%     (4) End of trial
% The number of samples within a trial is the maximum that is common to all
% trials. There is no per-trial temporal renormalization.
%
% Also provides the mouse's trajectory on each trial. Samples of the
% trajectory and traces are temporally aligned.
%

% Concatenate trial index information across all days, so that the set of
% frames common to all trials (over days) can be computed
trial_indices = [];
for di = md.valid_days
    trial_indices = [trial_indices; md.day(di).trial_indices]; %#ok<AGROW>
end
num_trials = size(trial_indices,1);
alignment_frames = trial_indices(:, align_idx);
[pre_offset, post_offset] = compute_frame_offsets(trial_indices,...
    1:num_trials, alignment_frames);
alignment_axis = pre_offset:post_offset;
num_common_frames = length(alignment_axis);

num_cells = md.num_cells;

X = zeros(num_cells, num_common_frames, num_trials);
xs = cell(num_trials,1);
ys = cell(num_trials,1);

for k = 1:num_trials
    % day and neuron indices
    di = trial_map(k,1); % Day index
    ni = md.get_indices(di);
    
    ti = trial_map(k,2); % Trial index
    traces = md.day(di).trials(ti).traces;
    
    % For each trial, compute the samples to keep
    tfi = md.day(di).trial_indices(ti,:); % trial frame indices
    tfi = tfi - (tfi(1)-1);
    af = tfi(align_idx); % alignment frame
    pre_frame = af + pre_offset;
    post_frame = af + post_offset;
    sampled_frames = pre_frame:post_frame;
    
    X(:,:,k) = traces(ni, sampled_frames);
    
    if md.day(di).is_tracking_loaded
        centroids = md.day(di).trials(ti).centroids;
        xs{k} = centroids(sampled_frames, 1);
        ys{k} = centroids(sampled_frames, 2);
    end
end