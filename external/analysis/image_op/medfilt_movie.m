function medfilt_movie(movie_in, varargin)
% Apply median filtering from file ('movie_in') to file ('movie_out').
% Optional argument specifies the half-width of the median filter.
%
% If 'movie_out' is left as an empty string, then default name will be
% provided.
%
% The medfilt parameters will also be saved in the output HDF5 file under
% the '/Medfilt' directory
%
% Inputs:
%   movie_in:  Name of incoming HDF5 movie
%
% Variable input arguments
%   movie_out: Name of outgoing HDF5 movie
%   medfilt_halfwidth: Half-width of the median filter kernel
%   'crop_edges': add it as an argument for cropping off zeros at the edges
%   introduced by median filter
%
% Example usage:
%   norm_movie('c9m7d12.hdf5','',15);
%


% Default dataset name for the movie
movie_dataset = '/Data/Images';

% Defaults
crop_edges = 0;
movie_out = '';

% Grab the movie parameters
[movie_size, in_data_type] = get_dataset_info(movie_in, movie_dataset);
height = movie_size(1);
width = movie_size(2);
num_frames = movie_size(3);


medfilt_halfwidth = 1;
if ~isempty(varargin)
    for k = 1:length(varargin)
        switch varargin{k}
            case 'movie_out'
                movie_out = varargin{k+1};
            case 'medfilt_halfwidth'
                medfilt_halfwidth = varargin{k+1};
            case 'crop_edges'
                crop_edges = 1;
        end
    end
    
end

if isempty(movie_out)
    [~, name] = fileparts(movie_in);
    movie_out = sprintf('%s_med.hdf5', name);
    fprintf('medfilt_movie: Output movie will be saved as "%s"\n', movie_out);
end
    

medfilt_neighborhood = (1+2*medfilt_halfwidth)*[1 1];

% Prepare output movie
output_height = height - 2*medfilt_halfwidth*crop_edges;
output_width  = width  - 2*medfilt_halfwidth*crop_edges;
h5create(movie_out, movie_dataset,...
         [output_height, output_width, num_frames],...
         'Datatype', in_data_type);
     
copy_hdf5_params(movie_in, movie_out);

h5create(movie_out, '/Medfilt/HalfWidth', 1, 'Datatype', 'int16');
h5write(movie_out, '/Medfilt/HalfWidth', int16(medfilt_halfwidth));

% Apply medfilt
frame_chunk_size = 1000;
[frame_chunks, num_chunks] = make_frame_chunks(num_frames, frame_chunk_size);

for i = 1:num_chunks
    fprintf('%s: Applying medfilt on frames %d to %d (out of %d)...\n',...
        datestr(now), frame_chunks(i,1), frame_chunks(i,2), num_frames);
    
    chunk_start = frame_chunks(i,1);
    chunk_count = frame_chunks(i,2) - frame_chunks(i,1) + 1;
    
    movie_chunk = h5read(movie_in, movie_dataset,...
                         [1 1 chunk_start],...
                         [height width chunk_count]);
    
    movie_chunk_med = zeros(output_height, output_width, chunk_count, in_data_type);
    
    for frame_idx = 1:size(movie_chunk,3)
        frame = movie_chunk(:,:,frame_idx);
        frame_med = medfilt2(frame, medfilt_neighborhood);
        
        if crop_edges == 1
            movie_chunk_med(:,:,frame_idx) = frame_med((1+medfilt_halfwidth):(end-medfilt_halfwidth),...
                                                   (1+medfilt_halfwidth):(end-medfilt_halfwidth));
        else
            movie_chunk_med(:,:,frame_idx) = frame_med;
        end
        
    end
    
    h5write(movie_out, movie_dataset,...
        movie_chunk_med,...
        [1 1 chunk_start],...
        [output_height output_width chunk_count]);
end
fprintf('%s: Done!\n', datestr(now));