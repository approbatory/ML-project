function zscore_movie(movie_in, movie_out, varargin)
% Compute pixel-wise z-scored movie from HDF5 file ('movie_in') to 
% file ('movie_out').
%
% If 'movie_out' is left as an empty string, then default name will be
% provided.
%
% The z-scoring parameters will also be saved in the output HDF5 file 
% under the '/zsc' directory
%
% Inputs:
%   movie_in:  Name of incoming HDF5 movie
%   movie_out: Name of outgoing HDF5 movie
%

if isempty(movie_out)
    [~, name] = fileparts(movie_in);
    movie_out = sprintf('%s_zsc.hdf5', name);
    fprintf('zscore_movie: Output movie will be saved as "%s"\n', movie_out);
end

% Default dataset name for the movie
movie_dataset = '/Data/Images';

% Grab the movie parameters
[movie_size, ~] = get_dataset_info(movie_in, movie_dataset);
height = movie_size(1);
width = movie_size(2);
num_frames = movie_size(3);

% Begin DFF processing
%------------------------------------------------------------
frame_chunk_size = 2500;
[frame_chunks, num_chunks] = make_frame_chunks(num_frames, frame_chunk_size);

% Compute the mean and std projetions
A = zeros(height, width);
A2 = zeros(height, width);
for i = 1:num_chunks
    fprintf('%s: Reading frames %d to %d for STD image (out of %d)...\n',...
        datestr(now), frame_chunks(i,1), frame_chunks(i,2), num_frames);

    chunk_start = frame_chunks(i,1);
    chunk_count = frame_chunks(i,2) - frame_chunks(i,1) + 1;

    movie_chunk = h5read(movie_in, movie_dataset,...
                         [1 1 chunk_start],...
                         [height width chunk_count]);

    for k = 1:size(movie_chunk,3)
        A = A + movie_chunk(:,:,k);
        A2 = A2 + movie_chunk(:,:,k).^2;
    end
end
A = A / num_frames;
A2 = A2 / num_frames;
S = sqrt(A2-A.^2);

imagesc(S);
colormap gray;
axis image;
title(sprintf('%s (STD)', strrep(movie_in, '_', '\_')));

input('zscore_movie: Press enter to proceed >> ');

% Prepare output movie
%------------------------------------------------------------
h5create(movie_out, movie_dataset,...
         [height width num_frames],...
         'ChunkSize', [height width 1],...
         'Datatype', 'single');
     
copy_hdf5_params(movie_in, movie_out);     

h5create(movie_out, '/zsc/S', size(S), 'Datatype', 'single');
h5write(movie_out, '/zsc/S', S);
h5create(movie_out, '/zsc/A', size(A), 'Datatype', 'single');
h5write(movie_out, '/zsc/A', A);

% Apply DFF
for i = 1:num_chunks
    fprintf('%s: Z-scoring frames %d to %d (out of %d)...\n',...
        datestr(now), frame_chunks(i,1), frame_chunks(i,2), num_frames);
    
    chunk_start = frame_chunks(i,1);
    chunk_count = frame_chunks(i,2) - frame_chunks(i,1) + 1;
    
    movie_chunk = h5read(movie_in, movie_dataset,...
                         [1 1 chunk_start],...
                         [height width chunk_count]);

	movie_chunk_zsc = zeros(size(movie_chunk), 'single');
    
    for frame_idx = 1:size(movie_chunk,3)
        movie_chunk_zsc(:,:,frame_idx) = (movie_chunk(:,:,frame_idx)-A)./S;
    end
    
    h5write(movie_out, movie_dataset,...
            movie_chunk_zsc,...
            [1 1 chunk_start],...
            [height width chunk_count]);
end
fprintf('%s: Done!\n', datestr(now));