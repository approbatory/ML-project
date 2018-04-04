PARLOOPS =  4;%64;


algs = my_algs({'mvnb2', 'ecoclin'}, {'original', 'shuf'});
dayset = auto_dayset('open_field'); dayset = dayset{1};
for ix = 1:numel(dayset)
    disp(dayset(ix).label);
    [ds, X, ks, errf] = load_day(dayset(ix),...
        'ds', {'nocells'},...
        'data', {'filling', 'binary', 'openfield', true});
    %spoof_samp_gen = indep_spoof(X, ks); %maybe change back, want to see
    %if indep_spoof is somehow a flawed approach
    spoof_samp_gen = @() shuffle(X,ks);
    for i = 1:numel(algs)
        [train_err{ix,i}, test_err{ix,i}] = evaluate_alg(algs(i),...
            X, ks, 'eval_f', errf,...
            'train_frac', 0.7, 'par_loops', PARLOOPS);
        [spoof_train_err{ix,i}, spoof_test_err{ix,i}] = evaluate_alg(algs(i),...
            spoof_samp_gen, ks, 'eval_f', errf,...
            'train_frac', 0.7, 'par_loops', PARLOOPS, 'X_is_func', true);
    end
end

%%
figure;
subplot(1,2,1);
plotmat('error', train_err, test_err, algs, dayset, 3);
subplot(1,2,2);
plotmat('error - on spoof', spoof_train_err, spoof_test_err, algs, dayset, 3);
%subplot(2,2,3);
%plotmat('Test error', test_err, algs, dayset);
%subplot(2,2,4);
%plotmat('Test error - on spoof', spoof_test_err, algs, dayset);