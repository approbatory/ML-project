function M = load_prairie_tif(path_to_tif)
% Prairie saves frames of a series as separate TIF files, where each file
% contains exactly one image. This loader simply enumerates all TIF files
% in a given directory, and return the result as a single movie matrix.
%
% Usage:
%   M = load_prairie_tif(pwd);
%

files = dir(fullfile(path_to_tif, '*.tif'));
num_frames = length(files);

% Open one image to get the X, Y dimensions
A = imread(fullfile(path_to_tif, files(1).name));
M = zeros(size(A,1), size(A,2), num_frames, class(A));

for i = 1:num_frames
    M(:,:,i) = imread(fullfile(path_to_tif, files(i).name));
end