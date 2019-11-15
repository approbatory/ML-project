function binned_indices = bin_frame_indices(orig_indices, bin_factor)
% Bin the frame indices temporally by 'bin_factor', but do not bin frames
% across trial boundaries. See also `bin_movie_in_time`.
%
% Inputs:
%   frame_indices: [num_trials x 4] matrix where the i-th row indicates the
%       frame indices of trial i as [start open-gate close-gate end]
%
%   bin_factor: Integer indicating the binning factor
%
% Output:
%   binned_frame_indices: Same dimensions as frame_indices, but where the
%       indices have been converted to account for temporal binning

[num_trials, K] = size(orig_indices);

frames_per_trial = orig_indices(:,end) - orig_indices(:,1) + 1;
binned_frames_per_trial = floor(frames_per_trial/bin_factor);

% Explicitly remove dangling frames from the original frame indices
for trial_idx = 1:num_trials
    max_index = orig_indices(trial_idx,1) + ...
                bin_factor * binned_frames_per_trial(trial_idx) - 1;
    orig_indices(trial_idx,:) = min(orig_indices(trial_idx,:),...
                                    max_index*ones(1,K)); % Apply clamp
end

in_trial_offsets = diff(orig_indices,1,2);
in_trial_offsets = cumsum(in_trial_offsets,2);

% Generate the binned indices
binned_start_indices = cumsum([1; binned_frames_per_trial(1:end-1)]);
binned_remaining_indices = repmat(binned_start_indices, 1, K-1) +...
                           floor(in_trial_offsets/bin_factor);

binned_indices = [binned_start_indices binned_remaining_indices];

end