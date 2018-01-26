algs = my_algs({'mvnb2', 'ecoclin'}, {'original', 'shuf'});
dayset = my_daysets('c14m4');
for ix = 1:numel(dayset)
    disp(dayset(ix).label);
    [ds, X, ks, errf] = load_day(dayset(ix),...
        'ds', {'deprobe', 'nocells'},...
        'data', {'filling', 'box'});
    subset = subset_moving(ds);
    for i = 1:numel(algs)
        [train_err{ix,i}, test_err{ix,i},...
            sub_train_err{ix,i}, sub_test_err{ix,i},...
            moving_train_err{ix,i}, moving_test_err{ix,i}] =...
            evaluate_alg(algs(i), X, ks, 'eval_f', errf,...
            'subset', subset, 'train_frac', 0.7, 'par_loops', 64);
    end
end

%%
subplot(2,3,1);
plotmat('Train error', train_err, algs, dayset);
subplot(2,3,2);
plotmat('Train error, on moving', sub_train_err, algs, dayset);
subplot(2,3,3);
plotmat('Moving data train error', moving_train_err, algs, dayset);

subplot(2,3,4);
plotmat('Test error', test_err, algs, dayset);
subplot(2,3,5);
plotmat('Test error, on moving', sub_test_err, algs, dayset);
subplot(2,3,6);
plotmat('Moving data test error', moving_test_err, algs, dayset);
%%
function plotmat(titl, errs, algs, dayset)
errnbar(cellfun(@mean, errs), cellfun(@(x) std(x)/sqrt(length(x)), errs));
set(gca, 'XTickLabel', {dayset.label});
set(gca, 'XTick', 1:numel(algs));
xtickangle(5);
legend(algs.name);
ylim([0 inf]);
ylabel('bin distance');
title(titl);
end