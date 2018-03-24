function [label, short, ds] = ds_autolabel(directory, name)
[label, short, ds] = ds_autolabel_core(directory, name);
label = [name(end-2:end) ', ' label];
end

function [label, short, ds] = ds_autolabel_core(directory, name)
if ischar(directory)
    ds = quick_ds(fullfile(directory, name), 'deprobe', 'nocells');
else
    ds = directory;
end
if ds.num_trials == 1
    label = 'open field';
    return;
end
goal_N = strcmp({ds.trials.goal}, 'north');
start_E = strcmp({ds.trials.start}, 'east');
goal_R = goal_N == start_E;

[label, short] = label_chunk(goal_N, goal_R);
if ~isempty(label)
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
short = [short_before short_after];
end

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