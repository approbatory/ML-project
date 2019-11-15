function bin_movie_in_time(movie_in, movie_out, bin_factor, trial_indices, movie_dataset)
% Temporally bin the provided movie M by bin_factor from HDF5 file ('movie_in')
% to file ('movie_out'). Does not bin frames across trial boundaries. 
% Dangling frames will be excised. Furthermore, non-trial frames -- if 
% they exist in the original movie -- will also be excised.
%
% If 'movie_out' is left as an empty string, then default name will be
% provided.
%
% The binning parameters will also be saved in the output HDF5 file under
% the '/TimeBin' directory
%
% Inputs:
%   movie_in:  Name of incoming HDF5 movie
%   movie_out: Name of outgoing HDF5 movie
%   bin_factor: Number of frames in the original movie that will correspond
%       to one frame in the output movie (integer)
%   trial_indices: [num_trials x 4] matrix whose i-th row indicates the
%       [start open-gate close-gate end] indices of Trial i.
%
% Example usage:
%   bin_movie_in_time('c9m7d12.hdf5','',2,frame_indices);
%

if isempty(movie_out)
    [~, name] = fileparts(movie_in);
    movie_out = sprintf('%s_ti%d.hdf5', name, bin_factor);
    fprintf('bin_movie_in_time: Output movie will be saved as "%s"\n', movie_out);
end

% Default dataset name for the movie
if ~exist('movie_dataset', 'var')
    movie_dataset = '/Data/Images';
end

% Grab the movie parameters
[movie_size, ~] = get_dataset_info(movie_in, movie_dataset);
height = movie_size(1);
width = movie_size(2);
num_frames = movie_size(3);

if (nargin < 4) % trial_indices not provided
    trial_indices = make_frame_chunks(num_frames, bin_factor * 1000);
end

assert(num_frames == trial_indices(end,end),...
       'Number of frames in movie does not match table of trial indices!');

% Begin temporal binning
%------------------------------------------------------------

% Get the binned indices
binned_indices = bin_frame_indices(trial_indices, bin_factor);
binned_frames_per_trial = binned_indices(:,end) - binned_indices(:,1) + 1;
num_binned_frames = sum(binned_frames_per_trial);
num_trials = size(binned_frames_per_trial,1);

% Prepare output movie
%------------------------------------------------------------
h5create(movie_out, movie_dataset,...
         [height width num_binned_frames],...
         'ChunkSize', [height width 1],...
         'Datatype', 'single');
     
copy_hdf5_params(movie_in, movie_out);

% Update the FPS
try
    fps = h5read(movie_in, '/Params/FrameRate');
    h5write(movie_out, '/Params/FrameRate', fps / bin_factor);
catch
    fprintf('Warning: Frame rate of file "%s" is unknown\n', movie_in);
end

h5create(movie_out, '/TimeBin/BinFactor', 1, 'Datatype', 'double');
h5write(movie_out, '/TimeBin/BinFactor', bin_factor);

for trial_idx = 1:num_trials
    fprintf('%s: Temporally binning Trial %d of %d...\n',...
        datestr(now), trial_idx, num_trials);
    
    % Retrieve frames from the original movie that will be used for the
    % binned output. Dangling frames are omitted!
    num_trial_frames = bin_factor * binned_frames_per_trial(trial_idx);
    
    trial_chunk = h5read(movie_in, movie_dataset,...
                         [1 1 trial_indices(trial_idx,1)],...
                         [height width num_trial_frames]);
    
    trial_chunk = reshape(trial_chunk, height, width,...
                          bin_factor, binned_frames_per_trial(trial_idx));
    trial_chunk = single(trial_chunk);
    
    trial_chunk = squeeze(mean(trial_chunk,3));
    
    h5write(movie_out, movie_dataset,...
            trial_chunk,...
            [1 1 binned_indices(trial_idx,1)],...
            [height width binned_frames_per_trial(trial_idx)]);
end
fprintf('%s: Done!\n', datestr(now));
