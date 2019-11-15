function trim_behavior_video(plusmaze_source, behavior_source, lick_source, trim, varargin)
% Extract the trial frames of the behavior video, so that the resulting
%   video lines up with the concatenated Miniscope movie (e.g. as produced
%   by `concatenateHDF5` or `concatenate_bigtiff`). The `trim` parameter
%   needs to match the values used in the Miniscope concatenation.
%
% Inputs:
%   plusmaze_source: Text file output from the plus maze
%   behavior_source: Behavior video (MPEG-4)
%   lick_source: Lickometer text file from the plus maze. Optional. Provide
%       empty matrix ('[]') to skip.
%   trim: Number of frames to drop from the beginning (trim[1]) and
%         end (trim[2]) of each trial
%   
% Optional input:
%   dropped_table: [N x 2] matrix where each row indicates the number of
%       dropped frames in a specified trial of the behavior video.
%
%       For example,
%    
%           dropped_table = [1 7; 22 6; 42 6];
%
%       indicates that Trial 1 has 7 dropped frames; Trial 22 has 6 dropped
%       frames; Trial 42 has 6 dropped frames.
%
% Example usage:
%   dropped_table = load('mouse7_day07_ego-left_dropped_behavior_frames.txt');
%   trim_behavior_video('mouse7_d07_ego-left.txt', 'mouse7_day07_ego-left.m4v', [15 5], dropped_table)
%



% Get the trial frames (beginning and end) according to PlusMaze output.
%   We regard the PlusMaze trial indices as ground truth!
%----------------------------------------------------------------------
orig_frame_indices = get_trial_frame_indices(plusmaze_source);
orig_frame_indices = orig_frame_indices(:,[1 4]); % Keep [Start End]

num_trials = size(orig_frame_indices,1);
num_frames = orig_frame_indices(num_trials,2); % Very last frame
fprintf('  PlusMaze output (%s) has %d frames\n', plusmaze_source, num_frames);

% Load the lickometer data. Note that licometer sampling is synchronized to
% the FPGA counter (and NOT the behavior video). Hence, we don't need to
% apply dropped behavior frame compensation.
%----------------------------------------------------------------------
is_lick_loaded = ~isempty(lick_source);
if is_lick_loaded
    lick_series = load(lick_source);
    num_lick_samples = length(lick_series);
    fprintf('  Lickometer series (%s) has %d samples\n', lick_source, num_lick_samples);
    assert(num_lick_samples == num_frames,...
           '  Number of samples (%d) in lickometer file is inconsistent with PlusMaze output (%d)!',...
           num_lick_samples, num_frames);
end

% Generate frame indices into the behavior video, with its dropped frames
%   (Optional parameter specifies a table of dropped frames)
%----------------------------------------------------------------------
dropped_frames = zeros(num_trials,1);
if ~isempty(varargin)
    dropped_frames_sparse = varargin{1}; % [Trial-number Num-dropped-frames]
    dropped_frames(dropped_frames_sparse(:,1)) = dropped_frames_sparse(:,2);
end

dropped_frame_indices = compute_dropped_frame_indices(orig_frame_indices, dropped_frames);

% Read behavior video, and make sure that the dropped trial indices
%   match up with the behavior video's number of frames
%----------------------------------------------------------------------
behavior_video = VideoReader(behavior_source);
num_behavior_frames = behavior_video.NumberOfFrames;
fprintf('  Behavior video (%s) has %d frames\n', behavior_source, num_behavior_frames);

assert(num_behavior_frames == dropped_frame_indices(end,2),...
       '  Dropped frame table is inconsistent with actual number of frames in behavior movie!');

% We don't handle the case where the number of frames in the behavior video
%   exceeds the frame count indicated by the PlusMaze
num_actual_dropped_frames = num_frames - num_behavior_frames;
assert(num_actual_dropped_frames >= 0,...
       '  Behavior video has more frames (%d) than indicated by PlusMaze output (%d)!',...
       num_behavior_frames, num_frames);

fprintf('  Behavior video shows %d dropped frames relative to PlusMaze output\n', num_actual_dropped_frames);

% Frames to keep
%----------------------------------------------------------------------
frames_to_keep = [orig_frame_indices(:,1)+trim(1) orig_frame_indices(:,2)-trim(2)];
num_trimmed_frames = sum(diff(frames_to_keep,1,2)+1); % Used for tracking progress

% Prepare output file
%----------------------------------------------------------------------
[~, name] = fileparts(behavior_source);
output_name = sprintf('%s_trim%d-%d', name, trim(1), trim(2));
trimmed_behavior_video = VideoWriter(output_name, 'MPEG-4');
trimmed_behavior_video.Quality = 100;
trimmed_behavior_video.FrameRate = 20; % FIXME: Don't hardcode
open(trimmed_behavior_video);

% Parametrize the "lick indicator"
lick_square_size = floor(behavior_video.Width / 10);
lick_square_border = 1;

write_idx = 0; % For tracking progress
for trial_idx = 1:num_trials
    frame_indices = frames_to_keep(trial_idx,1):frames_to_keep(trial_idx,2);
    
    % "Smear" the dropped frames within the trial
    behavior_frame_indices = interp1(orig_frame_indices(trial_idx,:),...
                                     dropped_frame_indices(trial_idx,:),...
                                     frame_indices,...
                                     'linear');
    behavior_frame_indices = round(behavior_frame_indices);
    
    for behavior_frame_idx = behavior_frame_indices
        A = read(behavior_video, behavior_frame_idx);
        A = rgb2gray(A);
        
        % Lick indicator at top right corner of video
        if is_lick_loaded
            lick_indicator = 255*ones(lick_square_size); % Max uint8
            if ~lick_series(behavior_frame_idx) % If no lick, mask white square with black
                lick_indicator((1+lick_square_border):(end-lick_square_border),...
                               (1+lick_square_border):(end-lick_square_border)) =...
                                    zeros(lick_square_size-2*lick_square_border);
            end
            A(1:lick_square_size, (end-(lick_square_size-1)):end) = lick_indicator;
        end
        
        writeVideo(trimmed_behavior_video, A);

        write_idx = write_idx + 1;
        if (mod(write_idx,1000)==0)
            fprintf('%s: Frames %d of %d written\n',...
                datestr(now), write_idx, num_trimmed_frames);
        end
    end
end
close(trimmed_behavior_video);

end % trim_behavior_video

function dropped_frame_indices = compute_dropped_frame_indices(orig_frame_indices, dropped_frames)
    % Recompute trial frame indices for indexing into the behavior video,
    %   taking into account the latter's dropped frames per each trial
    orig_trial_offsets = orig_frame_indices(2:end,1) - orig_frame_indices(1:(end-1),1);
    orig_trial_lengths = orig_frame_indices(:,2) - orig_frame_indices(:,1) + 1;

    dropped_trial_offsets = orig_trial_offsets - dropped_frames(1:(end-1));
    dropped_trial_lengths = orig_trial_lengths - dropped_frames;

    dropped_trial_starts = cumsum([orig_frame_indices(1,1); dropped_trial_offsets]);
    dropped_trial_ends = dropped_trial_starts + dropped_trial_lengths - 1;
    dropped_frame_indices = [dropped_trial_starts dropped_trial_ends];
end