PAR_LOOPS = 2;
algs = my_algs({'mvnb2', 'ecoclin'}, {'original', 'shuf'});
dayset = my_daysets('c11m1');
for ix = 1:numel(dayset)
    disp(dayset(ix).label);
    [ds, X, ks, errf] = load_day(dayset(ix),...
        'ds', {'deprobe', 'nocells'},...
        'data', {'filling', 'binary'});
    %spoof_samp_gen = indep_spoof(X, ks);
    spoof_samp_gen = @() shuffle(X,ks);
    subset = subset_moving(ds);
    for i = 1:numel(algs)
        [train_err{ix,i}, test_err{ix,i},...
            sub_train_err{ix,i}, sub_test_err{ix,i},...
            moving_train_err{ix,i}, moving_test_err{ix,i}] =...
            evaluate_alg(algs(i), X, ks, 'eval_f', errf,...
            'subset', subset, 'train_frac', 0.7, 'par_loops', PAR_LOOPS);
        [spoof_train_err{ix,i}, spoof_test_err{ix,i},...
            spoof_sub_train_err{ix,i}, spoof_sub_test_err{ix,i},...
            spoof_moving_train_err{ix,i}, spoof_moving_test_err{ix,i}] =...
            evaluate_alg(algs(i), spoof_samp_gen, ks, 'eval_f', errf,...
            'subset', subset, 'train_frac', 0.7, 'par_loops', PAR_LOOPS, 'X_is_func', true);
    end
end

%%
figure;
subplot(3,2,1);
plotmat('error', train_err, test_err, algs, dayset, 1.7);
subplot(3,2,3);
plotmat('error, on moving', sub_train_err, sub_test_err, algs, dayset, 1.3);
subplot(3,2,5);
plotmat('Moving data error', moving_train_err, moving_test_err, algs, dayset, 0.8);

subplot(3,2,2);
plotmat('Test error - on spoof', spoof_train_err, spoof_test_err, algs, dayset, 1.7);
subplot(3,2,4);
plotmat('Test error, on moving - on spoof', spoof_sub_train_err, spoof_sub_test_err, algs, dayset, 1.3);
subplot(3,2,6);
plotmat('Moving data test error - on spoof', spoof_moving_train_err, spoof_moving_test_err, algs, dayset, 0.8);