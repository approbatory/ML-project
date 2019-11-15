function M = load_movie_from_hdf5(source, varargin)
% Load a movie matrix [height x width x num_frames] from a HDF5 file.
%
% Optionally specify the frame range as a segment [start ... end]. Multiple
% segments are allowed as separate rows.
%
p = inputParser;
p.addOptional('segments', []);
p.addParameter('movie_dataset', '/Data/Images');
p.parse(varargin{:});
movie_dataset = p.Results.movie_dataset;

if ~isempty(p.Results.segments) % Read a subset of the movie
    segments = p.Results.segments; % [Start ... End]
    num_segments = size(segments,1);
    num_frames_per_segment = segments(:,end) - segments(:,1) + 1;
    
    [movie_size, movie_type] = get_dataset_info(source, movie_dataset);
    M = zeros(movie_size(1), movie_size(2), sum(num_frames_per_segment), movie_type); % preallocate
    
    % Index into movie matrix
    idx = 1;
    for k = 1:num_segments
        % Index into HDF5 file
        frame_start = segments(k,1);
        frame_count = num_frames_per_segment(k);

        inds = idx:(idx+num_frames_per_segment(k)-1);
        M(:,:,inds) = h5read(source, movie_dataset,...
                   [1 1 frame_start],...
                   [movie_size(1) movie_size(2) frame_count]);
        
        idx = idx + num_frames_per_segment(k);
    end
    
else % Read the whole file
    M = h5read(source, movie_dataset);
end