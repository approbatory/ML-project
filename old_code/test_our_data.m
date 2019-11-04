my_dayset = auto_dayset('open_field');
my_dayset = my_dayset{1};

algs = [my_algs({'gnb', 'lda'}),my_algs('ecoclin', {'shuf', 'original'})'];

for ix = 1:numel(my_dayset)
    %y = E_T{ix}.tracesEvents.position; %(T, 2)
    %if max(max(y)) > 1e3
    %    continue;
    %end
    %X = E_T{ix}.tracesEvents.rawTraces; %(T, C)
    ds = quick_ds(fullfile(my_dayset(ix).directory, my_dayset(ix).day), 'nocells');
    y = ds.trials.centroids;
    X = ds_dataset(ds, 'openfield', true, 'filling', 'traces', 'sparsify', false);
    
    for i = 1:numel(algs)
        [tr_err{ix,i}, te_err{ix,i}, ~] = evala(algs(i),...
            X, y, @(y) gen_place_bins(y,8,46),...
            'split', 'nonlocal', 'repeats', 4, 'verbose', true, 'use_par', true);
        [tr_err_sh{ix,i}, te_err_sh{ix,i}, ~] = evala(algs(i),...
            X, y, @(y) gen_place_bins(y,8,46),...
            'split', 'nonlocal', 'repeats', 4, 'verbose', true, 'shufboth', true, 'use_par', true);
        fprintf('Done ix=%d i=%d, te_err=%.2f +- %.2fcm\n', ix,i,mean(te_err{ix,i}), std(te_err{ix,i})/sqrt(length(te_err{ix,i})));
    end
end

%%
%up_to = numel(my_dayset);
figure;
berr({my_dayset.meta}, te_err(1:up_to,:), tr_err(1:up_to,:), {algs.short});
xlabel('Open field sessions');
ylabel('RMS error (cm)');
title('Place decoding error using traces (nonlocal split)');
ylim([0 25]);
figure;
berr({my_dayset.meta}, te_err_sh(1:up_to,:), tr_err_sh(1:up_to,:), {algs.short});
xlabel('Open field sessions');
ylabel('RMS error (cm)');
title('Place decoding error, trained/tested on shuffle');
ylim([0 25]);

%%
figure;
bsep({'mPFC open field'}, te_err(1,:), tr_err(1,:), 20, {{algs.short}}, 'RMS error (cm)');
%%
figure;
bsep({'mPFC open field, train/test shuffled'}, te_err_sh(1,:), tr_err_sh(1,:), 20, {{algs.short}}, 'RMS error (cm)');