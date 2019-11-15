function crop_movie(movie_in, movie_out, varargin)
% Crop HDF5 movie from file ('movie_in') to file ('movie_out'). Optional
% arguments can specify the cropping region programmatically. For example:
%   - crop_movie(movie_in, movie_out, 'trim', 5);
% will crop out 5 pixels from the borders of the movie.
%
% If no optional arguments are provided, then the cropping region is
% defined interactively.
%
% If 'movie_out' is left as an empty string, then default name will be
% provided.
%
% The cropping parameters will also be saved in the output HDF5 file under
% the '/Crop' directory
%
% Inputs:
%   movie_in:  Name of incoming HDF5 movie
%   movie_out: Name of outgoing HDF5 movie
%
% Example usage:
%   crop_movie('c9m7d12.hdf5','');
%

if isempty(movie_out)
    [~, name] = fileparts(movie_in);
    movie_out = sprintf('%s_cr.hdf5', name);
    fprintf('crop_movie: Output movie will be saved as "%s"\n', movie_out);
end

% Default dataset name for the movie
movie_dataset = '/Data/Images';
% Optional specification of the cropping ROI
x_bounds = [];
y_bounds = [];
use_projection = 0;
projection_type = 'min';

use_automin = 0;


if ~isempty(varargin)
    for k = 1:length(varargin)
        if ischar(varargin{k})
            vararg = lower(varargin{k});
            switch vararg
                case 'trim' % Trim borders
                    trim = varargin{k+1};
                    x_bounds = [1+trim width-trim];
                    y_bounds = [1+trim height-trim];
                case 'bounds' % Specify [x_bounds y_bounds]
                    bounds = varargin{k+1};
                    x_bounds = bounds(1:2);
                    y_bounds = bounds(3:4);
                case {'min', 'max'} % Compute projection image
                    use_projection = 1;
                    projection_type = vararg;
                case 'automin'
                    % Automatically computes the crop bounds, assuming 
                    % TurboReg registration (fills borders with 0's) with
                    % rotation _disabled_
                    use_automin = 1;
                    projection_type = 'min';
		case 'movie_dataset'
		    movie_dataset = varargin{k+1};
	        case 'chunk_size'
		    frame_chunk_size = varargin{k+1};
            end
        end
    end
end

% Grab the movie parameters
[movie_size, in_data_type] = get_dataset_info(movie_in, movie_dataset);
height = movie_size(1);
width = movie_size(2);
num_frames = movie_size(3);

% Begin crop processing
%------------------------------------------------------------

if ~exist('frame_chunk_size', 'var')
    frame_chunk_size = 2500;
end
[frame_chunks, num_chunks] = make_frame_chunks(num_frames, frame_chunk_size);

if use_projection || use_automin
    ref_frame = zeros(height, width);
    for i = 1:num_chunks
        fprintf('%s: Reading frames %d to %d for %s projection (out of %d)...\n',...
            datestr(now), frame_chunks(i,1), frame_chunks(i,2), projection_type, num_frames);
        
        chunk_start = frame_chunks(i,1);
        chunk_count = frame_chunks(i,2) - frame_chunks(i,1) + 1;
        
        movie_chunk = h5read(movie_in, movie_dataset,...
                             [1 1 chunk_start],...
                             [height width chunk_count]);
        
        % Actual operation depends on projection type
        switch projection_type
            case 'min'
                movie_chunk_min = min(movie_chunk, [], 3);
                if (i == 1) % Need to initialize
                    ref_frame = movie_chunk_min;
                else
                    ref_frame = min(ref_frame, movie_chunk_min);
                end

            case 'max'
                movie_chunk_max = max(movie_chunk, [], 3);
                if (i == 1)
                    ref_frame = movie_chunk_max;
                else
                    ref_frame = max(ref_frame, movie_chunk_max);
                end
        end
        
    end
    
    if use_automin
        [x_bounds, y_bounds] = compute_automin_bounds(ref_frame);
    end
else % Just display the middle frame
    ref_idx = floor(num_frames/2);
    ref_frame = h5read(movie_in, movie_dataset, [1 1 ref_idx], [height width 1]);
end

if isempty(x_bounds) % Cropping region was not specified programmatically
    imagesc(ref_frame);
    colormap gray;
    axis image;
    title(strrep(movie_in, '_', '\_'));
    fprintf('crop_movie: Please provide a rectangular region over the image.\n');
    fprintf('  Double click on the rectangle when done.\n');
    h_rect = imrect;
    rect_params = round(wait(h_rect));
    x_bounds = [rect_params(1) rect_params(1)+rect_params(3)-1];
    y_bounds = [rect_params(2) rect_params(2)+rect_params(4)-1];
end

% Show the cropped image
ref_frame_cropped = ref_frame(y_bounds(1):y_bounds(2),...
                              x_bounds(1):x_bounds(2));
imagesc(ref_frame_cropped);
colormap gray;
axis image;
title(sprintf('%s (cropped)',strrep(movie_in, '_', '\_')));

fprintf('crop_movie: Crop bounds are [x0 x1 y0 y1] = [%d %d %d %d]\n',...
    x_bounds(1), x_bounds(2), y_bounds(1), y_bounds(2));
input('crop_movie: Press enter to proceed >> ');

% Prepare output movie
%------------------------------------------------------------
cropped_height = y_bounds(2)-y_bounds(1)+1;
cropped_width = x_bounds(2)-x_bounds(1)+1;
h5create(movie_out, movie_dataset,...
         [cropped_height cropped_width, num_frames],...
         'ChunkSize', [cropped_height, cropped_width 1],...
         'Datatype', in_data_type); % Preserve data type

copy_hdf5_params(movie_in, movie_out);

h5create(movie_out, '/Crop/XBounds', [1 2], 'Datatype', 'double');
h5write(movie_out, '/Crop/XBounds', x_bounds);
h5create(movie_out, '/Crop/YBounds', [1 2], 'Datatype', 'double');
h5write(movie_out, '/Crop/YBounds', y_bounds);

for i = 1:num_chunks
    fprintf('%s: Cropping frames %d to %d (out of %d)...\n',...
        datestr(now), frame_chunks(i,1), frame_chunks(i,2), num_frames);
    
    chunk_start = frame_chunks(i,1);
    chunk_count = frame_chunks(i,2) - frame_chunks(i,1) + 1;
    
    movie_chunk = h5read(movie_in, movie_dataset,...
                         [1 1 chunk_start],...
                         [height width chunk_count]);
                         
    h5write(movie_out, movie_dataset,...
            movie_chunk(y_bounds(1):y_bounds(2), x_bounds(1):x_bounds(2), :),...
            [1 1 chunk_start],...
            [cropped_height cropped_width chunk_count]);
end
fprintf('%s: Done!\n', datestr(now));

end % crop_movie

function [x_bounds, y_bounds] = compute_automin_bounds(ref_frame)
    % Automatically compute the crop bounds, assuming TurboReg registration
    % (which fills borders with 0's) with rotation disabled.
    %
    % TODO: Generalize to account for rotation-enabled registration.
    
    [height, width] = size(ref_frame);
    
    % Coordinates for (roughly) the center of image
    x_mid = round(width/2);
    y_mid = round(height/2);
    
    % 1D slices through the image, then binarized
    x_cut = (ref_frame(y_mid, :) ~= 0);
    y_cut = (ref_frame(:, x_mid) ~= 0);
    
    x_bounds = [find(x_cut, 1, 'first') find(x_cut, 1, 'last')];
    y_bounds = [find(y_cut, 1, 'first') find(y_cut, 1, 'last')];
end
