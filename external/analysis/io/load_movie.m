function M = load_movie(source)
% Call the appropriate movie loading method based on the file extension
%   of the specified `source`
%
% 2015 01 28 Tony Hyun Kim

[~, ~, ext] = fileparts(source);
ext = lower(ext(2:end)); % Remove the leading dot

switch ext
    case {'tif', 'tiff'}
        M = load_movie_from_tif(source);
    case {'h5', 'hdf5'}
        M = load_movie_from_hdf5(source);
    otherwise
        % TODO: Handle error!
end