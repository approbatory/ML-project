function overwrite_frame_indices(plusmaze_source, new_frame_indices)
% Generate a new PlusMaze text file (appends "_new" to the filename)
%   with a new set of frame indices
%
% Inputs:
%   plusmaze_source: Name of the original text file (trial information
%       other that frame indices will be taken from this file)
%   new_frame_indices: [num_trials x 4] matrix where the i-th row indicates
%       the [start open-gate close-gate end] frames of trial i.
%
% Example usage:
%   overwrite_frame_indices('mouse7_day10_allo-south.txt', frame_indices);

% Read in information from the original PlusMaze source
[~, location_info, time] = parse_plusmaze(plusmaze_source);
num_trials = length(time);

% Generate a new PlusMaze output
[path_to, name, ~] = fileparts(plusmaze_source);
output_name = sprintf('%s_new.txt', name);
output_name = fullfile(path_to, output_name);

fid = fopen(output_name, 'w');
for i = 1:num_trials
    fprintf(fid, '%s %s %s %.3f %d %d %d %d\n',...
        location_info{i,1}, location_info{i,2}, location_info{i,3},...
        time(i),...
        new_frame_indices(i,1), new_frame_indices(i,2),...
        new_frame_indices(i,3), new_frame_indices(i,4));
end
fclose(fid);