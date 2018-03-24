% ds_hpc = quick_ds('../cohort14_dual/c14m6/c14m6d10', 'deprobe', 'nocells', 'cm', 'hpc_cm01_fix');
% ds_prl = quick_ds('../cohort14_dual/c14m6/c14m6d10', 'deprobe', 'nocells', 'cm', 'prl_cm01_fix');
% 
% alg = my_algs('linsvm', 0.1); %0.1 for L1 regularization
% 
% [X, ks] = ds_dataset(ds_prl,...
%     'selection', 0.3,...
%     'filling', 'box',...
%     'trials', strcmp({ds_prl.trials.start}, 'west'),...
%     'target', {ds_prl.trials.end});
% [train_error_prl, test_error_prl, model, fitinf] =...
%     evaluate_alg(alg,...
%     X, strcmp(ks, 'north'),...
%     'retain_models', true,...
%     'retain_fitinfo', true,...
%     'train_frac', 1);
% [~, order] = sort(-abs(model.Beta));
% order = order(1:nnz(model.Beta));
% see = @(n) histogram2(X(:, n), strcmp(ks, 'north')', 10);


ds_prl = cohort11{1}(6).ds;

alg = my_algs('linsvm', 0.1); %0.1 for L1 regularization

[X, ks] = ds_dataset(ds_prl,...
    'selection', 0.1,...
    'filling', 'traces',...
    'trials', strcmp({ds_prl.trials.start}, cohort11{1}(6).changing),...
    'target', {ds_prl.trials.end});
[train_error_prl, test_error_prl, model, fitinf] =...
    evaluate_alg(alg,...
    X, strcmp(ks, 'north'),...
    'retain_models', true,...
    'retain_fitinfo', true,...
    'train_frac', 1);
[vals, order] = sort(-abs(model.Beta));
order = order(1:nnz(model.Beta));
vals = -vals(1:nnz(model.Beta));
see = @(n) histogram2(X(:, n), strcmp(ks, 'north')', 10);

%%
best6 = order(1:6);
figure;
for n_ix = 1:6
    n = best6(n_ix);
    events_north = X(strcmp(ks, 'north'), n);
    events_south = X(strcmp(ks, 'south'), n);
    min_v = min(X(:,n));
    max_v = max(X(:,n));
    subplot(2,3,n_ix);
    histogram(events_north, 10, 'Normalization', 'pdf', 'BinLimits', [min_v max_v]);
    hold on;
    histogram(events_south, 10, 'Normalization', 'pdf', 'BinLimits', [min_v max_v]);
    title(sprintf('activity for neuron %d, ranked %d', n, n_ix));
    xlabel('trace value, no event detection');
    ylabel('normalized counts');
    legend north south
end
suplabel('c11m1d17 predictors for end arm at 0.1 fraction of turn', 't');