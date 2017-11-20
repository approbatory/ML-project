%This is the main script for end-arm decoding:

clear;
directory = '../c14m4'; %replace with wherever your data is

labels{1} = 'd15, ego left';
days{1} = 'c14m4d15';

labels{2} = 'd16, ego left -> allo south';
days{2} = 'c14m4d16';

labels{3} = 'd17, allo south';
days{3} = 'c14m4d17';

if ~exist('figs', 'dir')
    mkdir('figs');
end

%figure;
for i = 1:3
    ds = quick_ds(fullfile(directory, days{i}), 'deprobe', 'nocells');
    [poss, err, err_map] = decode_end_nb(ds, 0.005, 0.4, ...
                                         @gen_all_X_at_pos_closest_shuffled_start_keeped);
    fprintf('trained %s\n', labels{i});
    plot(poss, err, '-x');
    hold on;
    %view_err(ds, poss, err, err_map, labels{i});%, 'save', 'figs', 'hide');
end

xlabel('arm position');
ylabel('Multinomial NB err');
title('Multinomial NB errors vs. arm pos SHUFFLED (START ARM KEEPED)');
legend(labels{:});