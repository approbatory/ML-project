function norm_movie(movie_in, movie_out, varargin)
% Normalize HDF5 movie from file ('movie_in') to file ('movie_out').
% Each frame is divided by a disk-filtered version of itself. Optional
% argument specifies the radius of the disk-filter.
%
% If 'movie_out' is left as an empty string, then default name will be
% provided.
%
% The normalization parameters will also be saved in the output HDF5 file
% under the '/Norm' directory
%
% Inputs:
%   movie_in:  Name of incoming HDF5 movie
%   movie_out: Name of outgoing HDF5 movie
%
% Example usage:
%   norm_movie('c9m7d12.hdf5','',15);
%

if isempty(movie_out)
    [~, name] = fileparts(movie_in);
    movie_out = sprintf('%s_norm.hdf5', name);
    fprintf('norm_movie: Output movie will be saved as "%s"\n', movie_out);
end

% Default dataset name for the movie
movie_dataset = '/Data/Images';

% Grab the movie parameters
[movie_size, ~] = get_dataset_info(movie_in, movie_dataset);
height = movie_size(1);
width = movie_size(2);
num_frames = movie_size(3);

% Begin norm processing
%------------------------------------------------------------

% Optional parameters specify the disk filter radius
disk_radius = 15;
if ~isempty(varargin)
    disk_radius = varargin{1};
end
hDisk = fspecial('disk', disk_radius);

ref_idx = 1;
ref_frame = h5read(movie_in, movie_dataset, [1 1 ref_idx], [height width 1]);
ref_frame = imfilter(single(ref_frame), hDisk, 'replicate');
imagesc(ref_frame); axis image; colormap gray;
title(sprintf('%s (disk radius: %d)', strrep(movie_in,'_','\_'), disk_radius));

input('norm_movie: Press enter proceed >> ');

% Prepare output movie
%------------------------------------------------------------
h5create(movie_out, movie_dataset,...
         [height width num_frames],...
         'ChunkSize', [height width 1],...
         'Datatype', 'single');
     
copy_hdf5_params(movie_in, movie_out);     
     
h5create(movie_out, '/Norm/DiskRadius', 1, 'Datatype', 'double');
h5write(movie_out, '/Norm/DiskRadius', disk_radius);

% Apply norm
frame_chunk_size = 1000;
[frame_chunks, num_chunks] = make_frame_chunks(num_frames, frame_chunk_size);

for i = 1:num_chunks
    fprintf('%s: Normalizing %d to %d (out of %d)...\n',...
        datestr(now), frame_chunks(i,1), frame_chunks(i,2), num_frames);
    
    chunk_start = frame_chunks(i,1);
    chunk_count = frame_chunks(i,2) - frame_chunks(i,1) + 1;
    
    movie_chunk = h5read(movie_in, movie_dataset,...
                         [1 1 chunk_start],...
                         [height width chunk_count]);

	movie_chunk_norm = zeros(size(movie_chunk), 'single');
    
    for frame_idx = 1:size(movie_chunk,3);
        frame = single(movie_chunk(:,:,frame_idx));
        frame_f = imfilter(frame, hDisk, 'replicate');
        
        movie_chunk_norm(:,:,frame_idx) = frame ./ frame_f;
    end
    
    h5write(movie_out, movie_dataset,...
            movie_chunk_norm,...
            [1 1 chunk_start],...
            [height width chunk_count]);
end
fprintf('%s: Done!\n', datestr(now));
