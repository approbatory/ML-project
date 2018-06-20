ds = quick_ds('../open_field/om2-0405', 'nocells');

%using signif_filt

y = ds.trials.centroids;
X = ds_dataset(ds, 'openfield', true, 'filling', 'traces', 'sparsify', false);
X_filt = X(:, signif_filt);

binner = @(y) gen_place_bins(y, 8, 46);

algs = [my_algs({'gnb', 'lda'}),my_algs('ecoclin', {'shuf', 'original'})'];

for i = 1:numel(algs)
    [tr_err_plain{i}, te_err_plain{i}, ~] =...
        evala(algs(i), X, y, binner,...
        'split', 'nonlocal', 'repeats', 16, 'verbose', true);
    [tr_err_filt{i}, te_err_filt{i}, ~] =...
        evala(algs(i), X_filt, y, binner,...
        'split', 'nonlocal', 'repeats', 16, 'verbose', true);

    [tr_err_shuf{i}, te_err_shuf{i}, ~] =...
        evala(algs(i), X, y, binner,...
        'split', 'nonlocal', 'repeats', 16, 'verbose', true, 'shufboth', true);
    [tr_err_shfilt{i}, te_err_shfilt{i}, ~] =...
        evala(algs(i), X_filt, y, binner,...
        'split', 'nonlocal', 'repeats', 16, 'verbose', true, 'shufboth', true);
end

%%
tot_tr = [tr_err_plain;tr_err_filt;tr_err_shuf;tr_err_shfilt];
tot_te = [te_err_plain;te_err_filt;te_err_shuf;te_err_shfilt];

my_labels = {'original', 'filtered', 'pre-shuffle', 'pre-shuffle & filtered'};
figure;
berr(my_labels, tot_te, tot_tr, {algs.short});