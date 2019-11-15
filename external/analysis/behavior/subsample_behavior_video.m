function binned_indices = subsample_behavior_video(behavior_source, bin_factor, trial_indices)
% Temporally subsample the provided behavior video (MPEG-4) by 'bin_factor',
% while taking into account trial boundaries. See also `bin_movie_in_time'.
%
% Inputs:
%   - behavior_source: Name of behavior video (MPEG-4)
%   - bin_factor: Number of frames in the original movie that will
%       correspond to one frame in the binned movie.
%   - trial_indices: [num_trials x 4] matrix, where the i-th row indicates
%       the frame indices of Trial i as [start open-gate close-gate end]
%
% Outputs:
%   - binned_indices: [num_trials x 4] matrix where the i-th row indicates
%       the frame indices with respect to the binned movie.
%

% Get the binned indices
binned_indices = bin_frame_indices(trial_indices, bin_factor);
binned_frames_per_trial = binned_indices(:,end) - binned_indices(:,1) + 1;
num_trials = size(binned_frames_per_trial, 1);

% Input source
behavior_video = VideoReader(behavior_source);

% Output source
[~, name] = fileparts(behavior_source);
output_name = sprintf('%s_ti%d', name, bin_factor);
output_video = VideoWriter(output_name, 'MPEG-4');
output_video.Quality = 100;
output_video.FrameRate = behavior_video.FrameRate / bin_factor;
open(output_video);

for trial_idx = 1:num_trials
    % Read in trial at a time
    movie = read(behavior_video,...
                 [trial_indices(trial_idx,1) trial_indices(trial_idx,end)]);

    % Subsample the movie
    sub_inds = 1:binned_frames_per_trial(trial_idx);
    sub_inds = 1 + bin_factor*(sub_inds-1);
    sub_movie = movie(:,:,:,sub_inds);
    
    writeVideo(output_video, sub_movie);

    fprintf('%s: Trial %d of %d subsampled\n',...
        datestr(now), trial_idx, num_trials);
end
close(output_video);
