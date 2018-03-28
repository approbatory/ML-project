function [start, goal] = get_trial_info(direc)
run(fullfile(direc, 'data_sources.m'));
fid = fopen(fullfile(direc,ans.maze));
S = textscan(fid, '%s%s%s%f%d%d%d%d');
fclose(fid);

start = S{1};
goal = S{2};

probe_filter = strcmp(goal, 'east') | strcmp(goal, 'west');
start = start(~probe_filter);
goal = goal(~probe_filter);
end