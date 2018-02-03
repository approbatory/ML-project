%ds = quick_ds('../c14m4/c14m4d16', 'deprobe', 'nocells'); % for hpc
ds = quick_ds('../c11m1/c11m1d13', 'deprobe', 'nocells'); % for mpfc


trials_chosen = strcmp({ds.trials.start}, 'west'); % using the changing path
indices_chosen = 1:150;
indices_chosen = indices_chosen(trials_chosen);
points = 0:0.1:0.5;
tr = cell(numel(algs), numel(points)); te = tr;
algs = my_algs({'mvnb2', 'linsvm'}, {'original', 'shuf'});
for j = 1:length(points)
    [X_, ks_, em_] = ds_dataset(ds, 'selection', points(j), 'filling', 'binary',...
        'trials', trials_chosen, 'target', {ds.trials.end});
    for i = 1:numel(algs)
        alg = algs(i);
        [tr{i,j}, te{i,j}] = evaluate_alg(alg, X_, strcmp(ks_,'north'),...
            'eval_f', @(k,p) mean(~(k(:)==p(:))), 'par_loops', 2048);
    end
end
%%

tr_mean = cellfun(@mean, tr);
te_mean = cellfun(@mean, te);

tr_un = cellfun(@(x) std(x)/sqrt(length(x)), tr);
te_un = cellfun(@(x) std(x)/sqrt(length(x)), te);

colors = 'rgbcymk';

figure;
for i = 1:numel(algs)
    errorbar(points, tr_mean(i,:), tr_un(i,:), ['-.' colors(i)], 'DisplayName', [algs(i).name ' train']);
    hold on;
    errorbar(points, te_mean(i,:), te_un(i,:), ['-' colors(i)], 'DisplayName', [algs(i).name ' test']);
    hold on;
end
legend
xlabel('Fraction of turn completed');
ylabel('End arm classification error');
title('End arm prediction at positions along turn');

%TODO: maybe apply on independent spoofed data?