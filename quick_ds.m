function ds = quick_ds(dirname, varargin)
nocells = false;
deprobe = false;
for k = 1:length(varargin)
    if ischar(varargin{k})
        switch varargin{k}
            case 'nocells'
                nocells = true;
            case 'deprobe'
                deprobe = true;
        end
    end
end

start_dir = pwd;
cd(dirname);
cm_name = 'cm01';
if ~exist(cm_name, 'dir')
    cm_name = 'cm01-fix';
end
sources = data_sources;
XY = load(sources.tracking);
fid = fopen(sources.maze);
S = textscan(fid, '%s%s%s%f%d%d%d%d');
fclose(fid);
%fields = {'start', 'goal', 'end',...
%    'time', 'trial_start', 'gate_open',...
%    'gate_close', 'trial_end'};
%disp(XY);

class_file_name = file_pattern(cm_name, 'class*.txt');
%class_file_name = ls(fullfile('cm01','class*.txt'));
%if strcmp(class_file_name(end),newline)
%    class_file_name = class_file_name(1:end-1); %removing trailing newline
%end
lines = strsplit(fileread(class_file_name), '\n');
goodcells = [];
for l=lines
    if ~isempty(l{1}) && isempty(strfind(l{1}, 'not'))
        goodcells = [goodcells, sscanf(l{1}, '%d')];
    end
end
ds.num_cells = length(goodcells);


traces_filters_filename = file_pattern(cm_name, 'rec*.mat');
%traces_filters_filename = ls(fullfile('cm01','rec*.mat'));
%if strcmp(traces_filters_filename(end), newline)
%    traces_filters_filename = traces_filters_filename(1:end-1);
%end
%traces_filters = load(traces_filters_filename);
%traces = traces_filters.traces(:,goodcells);
s = load(traces_filters_filename, 'traces');
traces = s.traces;
traces = traces(:,goodcells);
if ~nocells
    s = load(traces_filters_filename, 'filters');
    filters = s.filters;
    filters = filters(:,:,goodcells);
    clear traces_filters;
    for j = 1:ds.num_cells
        ds.cells(j,1) = struct('im', filters(:,:,j),...
            'mask', filters(:,:,j)~=0,...
            'com', image_cm(filters(:,:,j)), 'label', 'cell');
    end
end

events_fname = file_pattern(cm_name, 'events*.mat');
%events_fname = ls(fullfile('cm01','events*.mat'));
%if strcmp(events_fname(end), newline)
%    events_fname = events_fname(1:end-1);
%end
events = load(events_fname);
events.events = events.events(goodcells);
events = events.events;

ds.trial_indices = [S{5}, S{6}, S{7}, S{8}];
ds.num_trials = size(ds.trial_indices,1);
ds.full_num_frames = ds.trial_indices(end,end);
start_cell = S{1}; goal_cell = S{2}; end_cell = S{3}; time_arr = S{4};
for i = 1:ds.num_trials
    trial_begin_frame = ds.trial_indices(i,1);
    trial_end_frame = ds.trial_indices(i,end);
    
    trial_events = cell(ds.num_cells, 1);
    for j = 1:ds.num_cells
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
    ds.trials(i,1) = struct('start',start_cell{i},'goal',...
        goal_cell{i},'end',end_cell{i},'time',time_arr(i),...
        'correct', strcmp(goal_cell{i},end_cell{i}), 'turn', turn,...
        'centroids', XY(trial_begin_frame:trial_end_frame,:),...
        'traces', traces(trial_begin_frame:trial_end_frame,:)',...
        'events', {trial_events});
    %ds.trials(i,1).events = trial_events;
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
    mask(i) = strcmp(ds.trials(i).end, 'north') || strcmp(ds.trials(i).end, 'south');
end
ds_deprobed = ds;
ds_deprobed.trials = ds_deprobed.trials(mask);
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