function dff_movie(movie_in, movie_out, varargin)
% Compute DFF from HDF5 file ('movie_in') to file ('movie_out').
%
% If 'movie_out' is left as an empty string, then default name will be
% provided.
%
% The DFF parameters will also be saved in the output HDF5 file under the
% '/DFF' directory
%
% Inputs:
%   movie_in:  Name of incoming HDF5 movie
%   movie_out: Name of outgoing HDF5 movie
%
% Example usage:
%   dff_movie('c9m7d12_norm.hdf5','');
%
p = inputParser;
p.addParameter('movie_dataset', '/Data/Images', @ischar);
p.addParameter('frame_chunk_size', 2500, @isscalar);
p.addParameter('F0', []);
p.parse(varargin{:});

if isempty(movie_out)
    %[~, name] = fileparts(movie_in);
    movie_out = sprintf('%s_dff.hdf5', movie_in);
    fprintf('dff_movie: Output movie will be saved as "%s"\n', movie_out);
end

% Default dataset name for the movie
movie_dataset = p.Results.movie_dataset;

% Grab the movie parameters
[movie_size, movie_type] = get_dataset_info(movie_in, movie_dataset);
height = movie_size(1);
width = movie_size(2);
num_frames = movie_size(3);

% Begin DFF processing
%------------------------------------------------------------
frame_chunk_size = p.Results.frame_chunk_size;
[frame_chunks, num_chunks] = make_frame_chunks(num_frames, frame_chunk_size);

if ~isempty(p.Results.F0)
    F0_image = p.Results.F0;
else
    % Construct the F0 image
    F0_image = zeros(height, width, 'single');
    for i = 1:num_chunks
        fprintf('%s: Reading frames %d to %d for F0 (out of %d)...\n',...
            datestr(now), frame_chunks(i,1), frame_chunks(i,2), num_frames);

        chunk_start = frame_chunks(i,1);
        chunk_count = frame_chunks(i,2) - frame_chunks(i,1) + 1;

        movie_chunk = single(h5read(movie_in, movie_dataset,...
                             [1 1 chunk_start],...
                             [height width chunk_count]));

        movie_chunk_sum = sum(movie_chunk,3);
        F0_image = F0_image + movie_chunk_sum;
    end
    F0_image = single(F0_image) / num_frames;
end

imagesc(F0_image);
colormap gray;
axis image;
title(sprintf('%s (F0)', strrep(movie_in, '_', '\_')));

input('dff_movie: Press enter to proceed >> ');

% Prepare output movie
%------------------------------------------------------------
h5create(movie_out, movie_dataset,...
         [height width num_frames],...
         'ChunkSize', [height width 1],...
         'Datatype', 'single');
     
copy_hdf5_params(movie_in, movie_out);     
     
h5create(movie_out, '/DFF/F0', size(F0_image), 'Datatype', 'single');
h5write(movie_out, '/DFF/F0', F0_image);

% Apply DFF
for i = 1:num_chunks
    fprintf('%s: Computing DFFs for frames %d to %d (out of %d)...\n',...
        datestr(now), frame_chunks(i,1), frame_chunks(i,2), num_frames);
    
    chunk_start = frame_chunks(i,1);
    chunk_count = frame_chunks(i,2) - frame_chunks(i,1) + 1;
    
    movie_chunk = h5read(movie_in, movie_dataset,...
                         [1 1 chunk_start],...
                         [height width chunk_count]);

	movie_chunk_dff = zeros(size(movie_chunk), 'single');
    
    for frame_idx = 1:size(movie_chunk,3)  
        movie_chunk_dff(:,:,frame_idx) = ...
            (single(movie_chunk(:,:,frame_idx))-F0_image) ./ F0_image;
    end
    
    h5write(movie_out, movie_dataset,...
            movie_chunk_dff,...
            [1 1 chunk_start],...
            [height width chunk_count]);
end
fprintf('%s: Done!\n', datestr(now));