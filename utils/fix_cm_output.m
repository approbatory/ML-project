origin = 'cm01';
fixed_output = 'cm01-fix';
rec_pattern = 'rec_*.mat';
class_pattern = 'class_*.txt';

if ~exist(origin, 'dir')
    error('origin dir: %s, does not exist', origin);
end
mkdir(fixed_output);

filt = strcmp({ds.cells.label}, 'cell');
rec_file = file_pattern(origin, rec_pattern);
rec_file_name = file_pattern(origin, rec_pattern, true);
S = load(rec_file);
S.traces = S.traces(:,filt);
S.filters = S.filters(:,:,filt);
S.info.num_pairs = sum(filt);
new_rec_fname = fullfile(fixed_output, rec_file_name);
save(new_rec_fname, '-struct', 'S');

n_cells = sum(filt);
class_file_name = file_pattern(origin, class_pattern, true);
new_class_file_name = fullfile(fixed_output, class_file_name);
fid = fopen(new_class_file_name, 'w');
for ix = 1:n_cells
    fprintf(fid, '%d, cell\n', ix);
end
fclose(fid);