function bin_movie_in_space(movie_in, movie_out, bin_factor, movie_dataset, frame_chunk_size)
% Temporally bin the provided movie M by bin_factor from HDF5 file ('movie_in')
% to file ('movie_out').
%
% If 'movie_out' is left as an empty string, then default name will be
% provided.
%
% The binning parameters will also be saved in the output HDF5 file under
% the '/SpaceBin' directory
%
% Inputs:
%   movie_in:  Name of incoming HDF5 movie
%   movie_out: Name of outgoing HDF5 movie
%   bin_factor: Pixel downsampling factor (integer!)
%

if isempty(movie_out)
    [~, name] = fileparts(movie_in);
    movie_out = sprintf('%s_sp%d.hdf5', name, bin_factor);
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

% Begin spatial binning
%------------------------------------------------------------
binned_height = floor(height / bin_factor);
binned_width = floor(width / bin_factor);

height_to_bin = bin_factor * binned_height;
width_to_bin = bin_factor*binned_width;

% Prepare output movie
%------------------------------------------------------------
h5create(movie_out, movie_dataset,...
         [binned_height binned_width num_frames],...
         'ChunkSize', [binned_height binned_width 1],...
         'Datatype', 'single');
     
copy_hdf5_params(movie_in, movie_out);

h5create(movie_out, '/SpaceBin/BinFactor', 1, 'Datatype', 'double');
h5write(movie_out, '/SpaceBin/BinFactor', bin_factor);

% Apply binning
if ~exist('frame_chunk_size', 'var')
    frame_chunk_size = 2500;
end
[frame_chunks, num_chunks] = make_frame_chunks(num_frames, frame_chunk_size);

for i = 1:num_chunks
    fprintf('%s: Spatially binning %d to %d (out of %d)...\n',...
        datestr(now), frame_chunks(i,1), frame_chunks(i,2), num_frames);
    
    chunk_start = frame_chunks(i,1);
    chunk_count = frame_chunks(i,2) - frame_chunks(i,1) + 1;
    
    movie_chunk = h5read(movie_in, movie_dataset,...
                         [1 1 chunk_start],...
                         [height width chunk_count]);

	movie_chunk = movie_chunk(1:height_to_bin, 1:width_to_bin, :);
    
    % Bin by reshaping
    movie_chunk = movie_chunk(:);
    movie_chunk = sum(reshape(movie_chunk, bin_factor, []), 1); % Bin along column
    movie_chunk = sum(reshape(movie_chunk, binned_height, bin_factor, []), 2); % Bin along rows
    movie_chunk = movie_chunk / bin_factor^2;
    
    movie_chunk = reshape(movie_chunk, binned_height, binned_width, chunk_count);

    h5write(movie_out, movie_dataset,...
            movie_chunk,...
            [1 1 chunk_start],...
            [binned_height binned_width chunk_count]);
end
fprintf('%s: Done!\n', datestr(now));
