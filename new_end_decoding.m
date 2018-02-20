%ds = quick_ds('../c14m4/c14m4d16', 'deprobe', 'nocells'); region = 'HPC'; % for hpc
ds = quick_ds('../c11m1/c11m1d13', 'deprobe', 'nocells'); region = 'mPFC';% for mpfc

PARLOOPS = 512;

trials_chosen = strcmp({ds.trials.start}, 'west'); % using the changing path (west for both)
indices_chosen = 1:150;
indices_chosen = indices_chosen(trials_chosen);
points = 0:0.1:0.5;
algs = my_algs({'mvnb2', 'linsvm'}, {'original', 'shuf'});
tr = cell(numel(algs), numel(points)); te = tr; tr_sh = tr; te_sh = tr; 
for j = 1:length(points)
    [X_, ks_, em_] = ds_dataset(ds, 'selection', points(j), 'filling', 'binary',...
        'trials', trials_chosen, 'target', {ds.trials.end});
    for i = 1:numel(algs)
        alg = algs(i);
        [tr{i,j}, te{i,j}] = evaluate_alg(alg, X_, strcmp(ks_,'north'),...
            'eval_f', @(k,p) mean(~(k(:)==p(:))), 'par_loops', PARLOOPS);
        [tr_sh{i,j}, te_sh{i,j}] = evaluate_alg(alg, @() shuffle(X_,strcmp(ks_,'north')), strcmp(ks_,'north'),...
            'eval_f', @(k,p) mean(~(k(:)==p(:))), 'X_is_func', true, 'par_loops', PARLOOPS);
    end
end
%%

tr_mean = cellfun(@mean, tr);
te_mean = cellfun(@mean, te);

tr_sh_mean = cellfun(@mean, tr_sh);
te_sh_mean = cellfun(@mean, te_sh);


tr_un = cellfun(@(x) std(x)/sqrt(length(x)), tr);
te_un = cellfun(@(x) std(x)/sqrt(length(x)), te);

tr_sh_un = cellfun(@(x) std(x)/sqrt(length(x)), tr_sh);
te_sh_un = cellfun(@(x) std(x)/sqrt(length(x)), te_sh);


colors = 'rgbcymk';

figure;
subplot(1,2,1);
for i = 1:numel(algs)
    errorbar(points, tr_mean(i,:), tr_un(i,:), ['-.' colors(i)], 'DisplayName', [algs(i).name ' train']);
    hold on;
    errorbar(points, te_mean(i,:), te_un(i,:), ['-' colors(i)], 'DisplayName', [algs(i).name ' test']);
    hold on;
end
ylim([0 0.15]);
legend
xlabel('Fraction of turn completed');
ylabel('End arm classification error');
title(sprintf('End arm prediction at positions along turn (%s)', region));

subplot(1,2,2);
for i = 1:numel(algs)
    errorbar(points, tr_sh_mean(i,:), tr_sh_un(i,:), ['-.' colors(i)], 'DisplayName', [algs(i).name ' train']);
    hold on;
    errorbar(points, te_sh_mean(i,:), te_sh_un(i,:), ['-' colors(i)], 'DisplayName', [algs(i).name ' test']);
    hold on;
end
ylim([0 0.15]);
legend
xlabel('Fraction of turn completed');
ylabel('End arm classification error');
title(sprintf('End arm prediction at positions along turn (shuffles of %s)', region));