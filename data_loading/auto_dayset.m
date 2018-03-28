function daysets = auto_dayset(names, varargin)


p = inputParser;
p.addRequired('names', @(x) iscell(x) || ischar(x));
p.addParameter('directory', '..', @ischar);
p.parse(names, varargin{:});
names = p.Results.names;

if ~iscell(names)
    names = {names};
end
daysets = cell(1,numel(names));
for ix = 1:numel(names)
    name = names{ix};
    directory = fullfile(p.Results.directory, name);
    S_files = dir(fullfile(directory,'c*m*d*'));
    if isempty(S_files)
        error('No c*m*d* data found under %s', directory);
    end
    for d_ix = 1:numel(S_files)
        s = S_files(d_ix);
        [label, short] = ds_autolabel(s.folder, s.name);
        daysets{ix}(d_ix) = struct('directory', s.folder,...
            'day', s.name, 'label', label, 'short', short,...
            'changing', get_changing(short));
    end
end

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
            label = [name(end-2:end) ', ' label];
            short = 'O';
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
        
        function [start, goal] = get_trial_info(direc)
            run(fullfile(direc, 'data_sources.m'));
            fid = fopen(fullfile(direc, ans.maze));
            S = textscan(fid, '%s%s%s%f%d%d%d%d');
            fclose(fid);
            
            start = S{1};
            goal = S{2};
            
            probe_filter = strcmp(goal, 'east') | strcmp(goal, 'west');
            start = start(~probe_filter);
            goal = goal(~probe_filter);
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
    end
    function ch = get_changing(short)
        switch short
            case {'RN', 'NR', 'LS', 'SL'}
                ch = 'west';
                return;
            case {'RS', 'SR', 'LN', 'NL'}
                ch = 'east';
                return;
        end
        ch = '';
    end
end