clear;

labels{1} = 'ego left (d15)'; days{1} = '~/brain/hpc/assets/c14m4d15/';
labels{2} = 'ego left -> allo south (d16)'; days{2} = '~/brain/hpc/assets/c14m4d16/';
labels{3} = 'allo south (d17)'; days{3} = '~/brain/hpc/assets/c14m4d17/';

for i = 1:3
    clear ds;
    ds = load_ds(days{i});
    [poss, err, err_map] = decode_end_nb(ds, 0.005, 0.4);
    view_err(ds, poss, err, err_map, labels{i}, 'save');
end