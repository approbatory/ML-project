MESSAGE = ['Running place decoding with shuffles, and shuffles on spoof, using'...
    ' GNB and ECOCSVM(linear, 1vsAll) on open field data, with traces filling using no regularization'...
    ' and using 4 train0.7/test0.3 splits, reporting average and sem'];



PARLOOPS =  4;%64;


%%algs = my_algs({'mvnb2', 'ecoclin_onevsall'}, {'original', 'shuf'}, true, 0);
algs = my_algs({'gnb', 'ecoclin_onevsall'}, {'original', 'shuf'}, true, 0);
dayset = auto_dayset('open_field'); dayset = dayset{1}([1 2 3 4 5 6 7 8 10 11 12 14]);
for ix = 1:numel(dayset)
    disp(dayset(ix).label);
    [ds, X, ks, errf] = load_day(dayset(ix),...
        'ds', {'nocells'},...
        'data', {'filling', 'traces', 'openfield', true, 'sparsify', false});    
        %'data', {'filling', 'copy_zeroed', 'openfield', true, 'sparsify', true});
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

save(sprintf('openfield_place_decoding_res_%s.mat', timestring),...
    'train_err', 'test_err', 'algs dayset', 'spoof_train_err',...
    'spoof_test_err', 'MESSAGE');
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