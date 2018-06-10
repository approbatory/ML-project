%comparing filling perf
PL = 16;


ds = quick_ds('../open_field/om3-0327', 'nocells');

fillings = {'binary', 'box', 'copy_zeroed', 'copy', 'traces'};

algs = my_algs({'ecoclin_onevsall'}, {'original', 'shuf'}, true, 0);


for ix = 1:numel(fillings)
    [X, ks, errf] = ds_dataset(ds, 'filling', fillings{ix},...
        'openfield', true, 'sparsify', ~strcmp(fillings{ix}, 'traces'));
    for i = 1:numel(algs)
        [train_err{ix,i}, test_err{ix,i}] = evaluate_alg(algs(i),...
            X, ks, 'eval_f', errf, 'train_frac', 0.7, 'par_loops', PL);
    end
    fprintf('Completed filling: %s\n', fillings{ix});
end

%%
fillings{3} = 'copy (zeroed)';
bin_cm = 5.12;
conv_to_cm = @(Q) cellfun(@(x) x.*5.12 ,Q, 'UniformOutput', false);

cm_train_err = conv_to_cm(train_err);
cm_test_err = conv_to_cm(test_err);
figure;
berr(fillings, cm_test_err, cm_train_err, {algs.name});

%plot labels:
title('Preprocessing Effect on Place Decoding');
xlabel('Trace preprocessing method');
ylabel('Mean bin distance, scaled by bin size (cm)');


dirname = 'graphs/place_decoding/open_field_hpc/filling_comparison';
fname = ['filling_comparison_' timestring '.png'];

set(gcf,'PaperUnits','inches','PaperPosition',2*[0 0 4 3]);
print(fullfile(dirname,fname), '-dpng', '-r100');