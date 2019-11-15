function save_movie_to_hdf5(movie, outname, varargin)
% Save a movie matrix [height x width x num_frames] into a HDF5 file.
%   Optional parameter to specify the dataset name.
%
% 2015 01 28 Tony Hyun Kim
if isempty(varargin)
    dataset_name = '/Data/Images';
else
    dataset_name = varargin{1};
end

[height, width, num_frames] = size(movie);
data_type = class(movie);

chunk_size = [height width 1];

h5create(outname, dataset_name, [height width num_frames],...
    'DataType', data_type,...
    'ChunkSize', chunk_size);
h5write(outname, dataset_name, movie);

% Create standard "/Params" directory
h5create(outname, '/Params/NumFrames', 1);
h5write(outname, '/Params/NumFrames', num_frames);