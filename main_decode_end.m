%This is the main script for end-arm decoding:

clear;
directory = '~/brain/hpc/assets'; %replace with wherever your data is

labels{1} = 'd15, ego left';
days{1} = 'c14m4d15';

labels{2} = 'd16, ego left -> allo south';
days{2} = 'c14m4d16';

labels{3} = 'd17, allo south';
days{3} = 'c14m4d17';

if ~exist('figs', 'dir')
    mkdir('figs');
end
for i = 1:3
    ds = quick_ds(fullfile(directory, days{i}), 'deprobe', 'nocells');
    [poss, err, err_map] = decode_end_nb(ds, 0.005, 0.4);
    view_err(ds, poss, err, err_map, labels{i}, 'save', 'figs', 'hide');
end