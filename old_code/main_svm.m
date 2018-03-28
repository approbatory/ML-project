clear;
rng(10);

directory = '../c14m4';

labels{1} = 'd15, ego left';
days{1} = 'c14m4d15';

labels{2} = 'd16, ego left to allo south';
days{2} = 'c14m4d16';

labels{3} = 'd17, allo south';
days{3} = 'c14m4d17';

if ~exist('figs', 'dir')
    mkdir('figs');
end

DO_SHUFFLE = true;
parfor i = 1:3
    ds = quick_ds(fullfile(directory, days{i}), 'deprobe', 'nocells');
    fprintf('loaded %s\n', labels{i});
    [poss{i}, err{i}, err_map{i}] = decode_end_svm(ds, 0.005, 0.4, DO_SHUFFLE);
    fprintf('trained %s\n', labels{i});
    %plot(poss, err, '-x');
    %hold on;
    %view_err(ds, poss{i}, err{i}, err_map{i}, labels{i}, 'SVM', 'save', 'newfigs', 'hide');
end

figure;
for i = 1:3
    plot(poss{i}, err{i}, '-x');
    hold on;
    %view_err(ds, poss{i}, err{i}, err_map{i}, labels{i}, 'save', 'figs', 'hide');
end

xlabel('arm position');
ylabel('SVM err');
plot_title = 'SVM errors vs. arm pos';
title_note = '';
if DO_SHUFFLE
    title_note = ' SHUFFLED';
end
title([plot_title title_note]);
legend(labels{:});
