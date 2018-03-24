function ds = quick_ds(dirname, varargin)
nocells = false;
deprobe = false;
cm_name_opt = '';
for k = 1:length(varargin)
    if ischar(varargin{k})
        switch varargin{k}
            case 'nocells'
                nocells = true;
            case 'deprobe'
                deprobe = true;
            case 'cm'
                cm_name_opt = varargin{k+1};
        end
    end
end

start_dir = pwd;
cd(dirname);

sources = data_sources;
XY = load(sources.tracking);
fid = fopen(sources.maze);
S = textscan(fid, '%s%s%s%f%d%d%d%d');
fclose(fid);



ds.trial_indices = [S{5}, S{6}, S{7}, S{8}];
ds.num_trials = size(ds.trial_indices,1);
ds.full_num_frames = ds.trial_indices(end,end);
ds.num_cells = [];

if ~isempty(cm_name_opt) && ~exist(cm_name_opt, 'dir')
    cd(start_dir);
    error('Directory for cellmax output not found: %s', cm_name_opt);
end
if exist(cm_name_opt, 'dir')
    [traces, events] = load_cm(cm_name_opt, nocells);
    ds.trials = get_trials(ds, traces, events, S, XY);
    ds.num_cells = size(ds.trials(1).traces,1);
elseif exist('hpc_cm01_fix', 'dir') && exist('prl_cm01_fix', 'dir')
    [traces_hpc, events_hpc] = load_cm('hpc_cm01_fix', nocells);
    trials_hpc = get_trials(ds, traces_hpc, events_hpc, S, XY);
    [traces_prl, events_prl] = load_cm('prl_cm01_fix', nocells);
    trials_prl = get_trials(ds, traces_prl, events_prl, S, XY);
    ds.num_cells = [size(trials_hpc(1).traces,1), size(trials_prl(1).traces,1)];
    ds.trials = [trials_hpc, trials_prl]; %%%%TODO CHECK!!!!
else
    cm_name = 'cm01-fix';
    if ~exist(cm_name, 'dir')
        cm_name = 'cm01';
    end
    [traces, events] = load_cm(cm_name, nocells);
    ds.trials = get_trials(ds, traces, events, S, XY);
    ds.num_cells = size(ds.trials(1).traces,1);
end

cd(start_dir);

if deprobe
    ds = de_probe(ds);
end
end



function cm = image_cm(im)
[W,H] = size(im);
[X, Y] = meshgrid(1:H,1:W);
tot_weight = sum(sum(im));
mean_x = sum(sum(X.*im)) ./ tot_weight;
mean_y = sum(sum(Y.*im)) ./ tot_weight;
cm = [mean_x; mean_y];
end


function [ds_deprobed] = de_probe(ds)
%DE_PROBE Remove probe trials from DaySummary object.
n_trials = length(ds.trials);
mask = false(1,n_trials);
for i = 1:n_trials
    mask(i) = strcmp(ds.trials(i,1).end, 'north') || strcmp(ds.trials(i,1).end, 'south');
end
ds_deprobed = ds;
ds_deprobed.trials = ds_deprobed.trials(mask,:);
ds_deprobed.trial_indices = ds_deprobed.trial_indices(mask,:);
ds_deprobed.num_trials = sum(mask);
end

function res = file_pattern(d, pat)
S = dir(fullfile(d,pat));
if numel(S) ~= 1
    error('There must be exactly one file matching the pattern: %s/%s', d, pat);
end
res = fullfile(S.folder, S.name);
end

function exists_one = was_found(d, pat)
S = dir(fullfile(d,pat));
exists_one = numel(S) == 1;
end

function [traces, events, varargout] = load_cm(cm_name, nocells)
filtering = was_found(cm_name, 'class*.txt');
if filtering
class_file_name = file_pattern(cm_name, 'class*.txt');
lines = strsplit(fileread(class_file_name), '\n');
goodcells = [];
for l=lines
    if ~isempty(l{1}) && isempty(strfind(l{1}, 'not'))
        goodcells = [goodcells, sscanf(l{1}, '%d')];
    end
end
num_cells = length(goodcells);
end
traces_filters_filename = file_pattern(cm_name, 'rec*.mat');
s = load(traces_filters_filename, 'traces');
traces = s.traces;
if filtering
traces = traces(:,goodcells);
else
    num_cells = size(traces,2);
    goodcells = true(1,num_cells);
end
if ~nocells
    s = load(traces_filters_filename, 'filters');
    filters = s.filters;
    filters = filters(:,:,goodcells);
    clear traces_filters;
    for j = 1:num_cells
        cells(j,1) = struct('im', filters(:,:,j),...
            'mask', filters(:,:,j)~=0,...
            'com', image_cm(filters(:,:,j)), 'label', 'cell');
    end
    varargout{1} = cells;
end
events_fname = file_pattern(cm_name, 'events*.mat');
events = load(events_fname);
events.events = events.events(goodcells);
events = events.events;
end

function trials = get_trials(ds, traces, events, S, XY)
start_cell = S{1}; goal_cell = S{2}; end_cell = S{3}; time_arr = S{4};
num_cells = size(traces, 2);
for i = 1:ds.num_trials
    trial_begin_frame = ds.trial_indices(i,1);
    trial_end_frame = ds.trial_indices(i,end);
    
    trial_events = cell(num_cells, 1);
    for j = 1:num_cells
        b = trial_begin_frame; e = trial_end_frame;
        if(isempty(events(j).auto))
            events(j).auto = zeros(1,3);
        end
        fst = events(j).auto(:,1);
        snd = events(j).auto(:,2);
        
        relevance_mask = (((fst >= b) & (fst <= e))| isinf(fst)) & ((snd >= b) & (snd <= e));
        
        relevant_events = events(j).auto(relevance_mask,:);
        relevant_events(:,1) = relevant_events(:,1) - double(b)+1;
        relevant_events(:,2) = relevant_events(:,2) - double(b)+1;
        trial_events{j} = relevant_events;
    end
    
    
    turn = 'left';
    if strcmp(start_cell{i},'east') && strcmp(end_cell{i},'north')
        turn = 'right';
    end
    if strcmp(start_cell{i},'west') && strcmp(end_cell{i},'south')
        turn = 'right';
    end
    if strcmp(start_cell{i},'north') && strcmp(end_cell{i},'west')
        turn = 'right';
    end
    if strcmp(start_cell{i},'south') && strcmp(end_cell{i},'east')
        turn = 'right';
    end
    trials(i,1) = struct('start',start_cell{i},'goal',...
        goal_cell{i},'end',end_cell{i},'time',time_arr(i),...
        'correct', strcmp(goal_cell{i},end_cell{i}), 'turn', turn,...
        'centroids', XY(trial_begin_frame:trial_end_frame,:),...
        'traces', traces(trial_begin_frame:trial_end_frame,:)',...
        'events', {trial_events});
end

end