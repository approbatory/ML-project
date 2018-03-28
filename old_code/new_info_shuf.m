PARLOOPS = 1;
alg = my_algs('ecoclin');
alg_shuf = my_algs('ecoclin', 'shuf');
dayset = my_daysets('c14m4');
for ix = 1:numel(dayset)
    [ds, X, ks, errf] = load_day(dayset(ix), 'ds', {'deprobe', 'nocells'},...
        'data', {'filling', 'binary'});
    [tr_base{ix}, te_base{ix}] = evaluate_alg(alg, X, ks, 'eval_f', errf, 'par_loops', PARLOOPS);
    report('Base', te_base{ix});
    [tr_shuf{ix}, te_shuf{ix}] = evaluate_alg(alg_shuf, X, ks, 'eval_f', errf, 'par_loops', PARLOOPS);
    report('Shuf', te_shuf{ix});
    rand_order = randperm(ds.num_cells);
    for n = 1:ds.num_cells
        fprintf('Cell %d of %d\n', n, ds.num_cells);
        shuf_mask = rand_order <= n;
        mod_alg = alg;
        mod_alg.train = @(X,ks) alg.train(shuffle(X,ks,'subset',shuf_mask),ks);
        [tr{ix,n}, te{ix,n}] = evaluate_alg(mod_alg, X, ks, 'eval_f', errf, 'par_loops', PARLOOPS);
        report('Partials', te{ix,n});
    end
end
save NEW_INFO_SHUF.mat

function report(l, te)
fprintf('%s Err: %f +- %f\n', l, mean(te), std(te)/sqrt(length(te)));
end