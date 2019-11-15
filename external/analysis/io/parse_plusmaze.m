function [frame_indices,location_info,time] = parse_plusmaze(source)

% Parse information from plusmaze textfile 
% 
% Example arguments:
% source = '/Volumes/COHORT9/cohort9-herrin224/mouse5/day20_ego-right/mouse5_day20_ego-right.txt';
%
% 2015-02-03 Jessica Maxey

fid = fopen(source);
maze_data = textscan(fid, '%s %s %s %f %d %d %d %d');
fclose(fid);

frame_indices = [maze_data{5} maze_data{6} maze_data{7} maze_data{8}];
location_info = [maze_data{1} maze_data{2} maze_data{3}];
time = maze_data{4};