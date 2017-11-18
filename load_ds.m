function ds = load_ds(dir, leave_probes)
if nargin == 1
    leave_probes = false;
end
start_dir = pwd;
cd(dir);
sources = data_sources;
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
    end
end
cd(start_dir);
if ~leave_probes
    de_probe(ds);
end
end
