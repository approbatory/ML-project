function binned_indices = bin_behavior_video(behavior_source, bin_factor, trial_indices)
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
num_binned_frames = sum(binned_frames_per_trial);
num_trials = size(binned_frames_per_trial, 1);

% Input source
behavior_video = VideoReader(behavior_source);
data_type = class(read(behavior_video, 1)); % e.g. 'uint8'

% Output source
[~, name] = fileparts(behavior_source);
output_name = sprintf('%s_ti%d', name, bin_factor);
output_video = VideoWriter(output_name, 'MPEG-4');
output_video.Quality = 100;
output_video.FrameRate = behavior_video.FrameRate / bin_factor;
open(output_video);

write_idx = 0;
for trial_idx = 1:num_trials
    for k = 1:binned_frames_per_trial(trial_idx)
        % Indices into the original movie
        frame_start = trial_indices(trial_idx,1) + bin_factor*(k-1);
        frame_end   = frame_start + (bin_factor-1);
        
        frames = read(behavior_video, [frame_start frame_end]);
        frames = squeeze(frames(:,:,1,:)); % Convert RGB to grayscale
        frames = single(frames); % To avoid integer arithmetic on 'mean'
        
        % Compute the mean, then recast to the original data type
        mean_frame = cast(mean(frames,3), data_type);
        
        writeVideo(output_video, mean_frame);

        write_idx = write_idx + 1;
        
        if (mod(write_idx,1000)==0)
            fprintf('%s: Frames %d of %d written\n',...
                datestr(now), write_idx, num_binned_frames);
        end
    end
end
close(output_video);
