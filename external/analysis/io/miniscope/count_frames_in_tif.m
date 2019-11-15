function [tif_frames, tif_filenames] = count_frames_in_tif(path_to_tif)
% Computes the number of frames contained in all TIF files of
%   the specified directory
%
% Example use: "count_frames_in_tif(pwd)"
%   in a directory that contains TIF files
%

tif_files = dir(fullfile(path_to_tif,'*.tif'));
num_files = length(tif_files);
fprintf('Found %d TIF files in "%s"\n', num_files, path_to_tif);

tif_frames = zeros(num_files, 1);
tif_filenames = cell(num_files, 1);

total_frames = 0;
for i = 1:num_files
    tif_filename = fullfile(path_to_tif,tif_files(i).name);
    tif_info = imfinfo(tif_filename);
    
    num_frames = length(tif_info);
    total_frames = total_frames + num_frames;
    
    tif_filenames{i} = tif_filename;
    tif_frames(i) = num_frames;
    fprintf('  %d: "%s" has %d frames\n',...
        i, tif_filename, num_frames);
end

fprintf('  Total frame count is %d\n', total_frames);