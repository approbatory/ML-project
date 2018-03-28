function [label, short] = ds_autolabel(directory, name)
if ischar(directory)
    [start, goal] = get_trial_info(fullfile(directory, name));
    num_trials = length(goal);
else
    ds = directory;
    start = {ds.trials.start};
    goal = {ds.trials.goal};
    num_trials = ds.num_trials;
end
if num_trials == 1
    label = 'open field';
    return;
end
goal_N = strcmp(goal, 'north');
start_E = strcmp(start, 'east');
goal_R = goal_N == start_E;

[label, short] = label_chunk(goal_N, goal_R);
if ~isempty(label)
    label = [name(end-2:end) ', ' label];
    return;
end

switch_at = 50;
goal_N_before = goal_N(1:switch_at);
goal_N_after = goal_N(switch_at+1:end);

goal_R_before = goal_R(1:switch_at);
goal_R_after = goal_R(switch_at+1:end);

[label_before, short_before] = label_chunk(goal_N_before, goal_R_before);
[label_after, short_after] = label_chunk(goal_N_after, goal_R_after);

if isempty(label_before) || isempty(label_after)
    error('could not determine label for %s/%s', directory, name);
end
label = [label_before ' to ' label_after];
label = [name(end-2:end) ', ' label];
short = [short_before short_after];

    function [label, short] = label_chunk(goal_N, goal_R)
        if all(goal_N)
            label = 'allo north';
            short = 'N';
            return;
        end
        
        if all(~goal_N)
            label = 'allo south';
            short = 'S';
            return;
        end
        
        if all(goal_R)
            label = 'ego right';
            short = 'R';
            return;
        end
        
        if all(~goal_R)
            label = 'ego left';
            short = 'L';
            return;
        end
        label = '';
        short = '';
    end
end