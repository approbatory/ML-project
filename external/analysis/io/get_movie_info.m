function [movie_size, movie_type] = get_movie_info(hdf5_source)
% Get the size and type (e.g. uint16 / single) of a movie stored in HDF5

dataset_name = '/Data/Images';
[movie_size, movie_type] = get_dataset_info(hdf5_source, dataset_name);