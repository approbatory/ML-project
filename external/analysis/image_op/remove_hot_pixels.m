function remove_hot_pixels(movie_in,movie_out,varargin)
% Remove hot pixels in the raw movie.
%
% Input arguments:
%   movie_in:  Name of incoming HDF5 movie
%   movie_out: Name of outgoing HDF5 movie 
%
% Optional input argument:
%   num_hot_pixels: Maximum estimated number of hot pixels. Default is 30.
%
% Example usage:
%   remove_hot_pixels('c9m7d12.hdf5','');
%

if ~isempty(varargin)
    num_hot_pixels = varargin{1};
    if ~isnumeric(num_hot_pixels)
        error('Maximum estimated number of hot pixels must be a numerical value');
    end
else
    num_hot_pixels = 30;
end
fprintf('%s: Loading movie into memory... \n',datestr(now));

M = load_movie(movie_in);
[height,width,num_frames] = size(M);

fprintf('%s: Removing hot pixels... \n',datestr(now));

min_proj = double(min(M,[],3));
median_min_proj = medfilt2(min_proj,[3,3]);

%Difference between the min projection before and after median filtering
%determines hot pixels
diff_proj = min_proj - median_min_proj;

% Remove corner pixels introduced due to median filter
for y = [1,height]
    for x = [1,width]
        diff_proj(y,x) = 0;
    end
end

% Threshold for calling a pixel "hot"
q = quantile(diff_proj(:),1-num_hot_pixels/(height*width));

idx = find(diff_proj>q);
[hot_y,hot_x] = ind2sub([height,width],idx);

% Median filtering only for the hot pixels
for i = 1:length(hot_y)
    x_neighbors = max(min(hot_x(i) + (-1:1),width),1);
    y_neighbors = max(min(hot_y(i) + (-1:1),height),1);
    M_neighbors = reshape(M(y_neighbors,x_neighbors,:),...
        length(x_neighbors)*length(y_neighbors),num_frames);
    M(hot_y(i),hot_x(i),:) = median(M_neighbors,1);
end

fprintf('%s: Writing the output movie to disk... \n',datestr(now));

if isempty(movie_out)
    [~, name] = fileparts(movie_in);
    movie_out = sprintf('%s_rm.hdf5', name);
    fprintf('Output movie will be saved as "%s"\n', movie_out);
end

% Default dataset name for the movie
movie_dataset = '/Data/Images';

% Prepare output movie
%------------------------------------------------------------
h5create(movie_out, movie_dataset,...
         [height width num_frames],...
         'ChunkSize', [height width 1],...
         'Datatype', 'uint16');
     
copy_hdf5_params(movie_in, movie_out);  

h5write(movie_out, movie_dataset,M);

fprintf('%s: All done! \n',datestr(now));