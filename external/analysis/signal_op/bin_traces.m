function binned_traces = bin_traces(traces, orig_indices, bin_factor)
% Temporally bin the provided `traces` by `bin_factor`. Does not bin frames
% across trial boundaries.
%
% Inputs:
%   traces: [num_frames x num_traces] matrix containing the traces to be
%           binned. (Binning occurs along the column direction.)
%   orig_indices: [num_trials x 4] matrix containing the (unbinned) trial
%           frame indices.
%   bin_factor: Number of frames in the original trace that will correspond
%           to one frame in the output trace.
%

% Get the binned indices
binned_indices = bin_frame_indices2(orig_indices, bin_factor);
binned_frames_per_trial = binned_indices(:,end) - binned_indices(:,1) + 1;
num_binned_frames = sum(binned_frames_per_trial);
num_trials = size(binned_frames_per_trial, 1);

% Initialize binned traces
num_traces = size(traces, 2);
binned_traces = zeros(num_binned_frames, num_traces);

for trial_idx = 1:num_trials
    % Retrieve frames from the original movie that will be used for the
    % binned output. Dangling frames are omitted!
    num_trial_frames = bin_factor * binned_frames_per_trial(trial_idx);
    
    trial_start = orig_indices(trial_idx, 1);
    trial_end   = trial_start + num_trial_frames - 1;
    
    trial_chunk = traces(trial_start:trial_end, :);
    
    trial_chunk = reshape(trial_chunk,...
        bin_factor, binned_frames_per_trial(trial_idx), num_traces);
    
    trial_chunk = squeeze(mean(trial_chunk,1));
    
    % Save to binned output
    binned_start = binned_indices(trial_idx, 1);
    binned_end   = binned_start + binned_frames_per_trial(trial_idx) - 1;
    binned_traces(binned_start:binned_end, :) = trial_chunk;
end