function ds = load_ds(dir, leave_probes)
% Load a DaySummary object from a directory, by default removes
% probe-trials, unless specified in the second argument
if nargin == 1
    leave_probes = false;
end
start_dir = pwd;
cd(dir);
sources = data_sources;
%video is unnecessary:
%sources = rmfield(sources, 'behavior');
if exist('cm01-fix', 'dir')
    ds = DaySummary(sources, 'cm01-fix');
elseif exist('cm01', 'dir')
    ds = DaySummary(sources, 'cm01');
    class_file_name = ls('cm01/class*.txt');
    class_file_name = class_file_name(1:end-1); %removing trailing newline
    lines = strsplit(fileread(class_file_name), '\n');
    goodcells = [];
    for l=lines
        if ~isempty(l{1}) && isempty(strfind(l{1}, 'not'))
            goodcells = [goodcells, sscanf(l{1}, '%d')];
        end
    end
    ds.num_cells = length(goodcells);
    ds.cells = ds.cells(goodcells);
    for i=1:ds.num_trials
        ds.trials(i).events = ds.trials(i).events(goodcells);
        ds.trials(i).traces = ds.trials(i).traces(goodcells,:);
    end
end
cd(start_dir);
if ~leave_probes
    de_probe(ds);
end
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

