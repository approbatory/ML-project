clear;
mv_pre = @(X) int64(X ~= 0);
mv_train = @(X_train, ks_train) fitcnb(X_train, ks_train, 'Distribution', 'mvmn',...
    'CategoricalPredictors', 'all');
mv_test = @(model, X_test) predict(model, X_test);

ecoc_pre = @(X) X;
ecoc_train = @(X_train, ks_train) fitcecoc(X_train, ks_train);
ecoc_test = @(model, X_test) predict(model, X_test);

%%%choose model here: need to implement shuffled ecoc, compare with it
%pre = mv_pre; train = mv_train; test = mv_test;
pre = ecoc_pre; train = ecoc_train; test = ecoc_test;
%%
ds = quick_ds(fullfile('../c14m4', 'c14m4d15'), 'deprobe', 'nocells');
reindexed_trial_indices = zeros(size(ds.trial_indices));
closing = 0;
for i = 1:size(ds.trial_indices,1)
    r = ds.trial_indices(i,:);
    reindexed_trial_indices(i,:) = r - r(1) + closing + 1;
    closing = reindexed_trial_indices(i,end);
end

eval_over = false(reindexed_trial_indices(end,end),1);
for r = reindexed_trial_indices'
    eval_over(r(2):r(3)) = true;
end


%%
X = cell2mat(gen_place_decoding_X(ds));
[ks, D] = bin_space(cell2mat(preprocess_xy(ds)));
K = size(D,1);
training_frac = 0.7;

dist_func = @(k,p) D(sub2ind(size(D), k, p));
err_func = @(k,p) k~=p;
mean_dist = @(k,p) mean(dist_func(k,p));
mean_err = @(k,p) mean(err_func(k,p));
keep_it = @(k,p) [k,p];

err_funcs = {mean_dist, mean_err, @confusionmat, keep_it};
func_names = {'dist', '0-1', 'conf', 'keep'};

[main, sub] = ...
    evaluate_decoder(pre, train, test, X, ks, eval_over,...
    err_funcs, func_names, training_frac);
sub_only = evaluate_decoder(pre, train, test,...
    X(eval_over,:), ks(eval_over), [],...
    err_funcs, func_names, training_frac);

fprintf('main:\ttrain dist: %f\ttest dist: %f\n', main(1).train_err, main(1).test_err);
fprintf('main:\ttrain err: %f\ttest err: %f\n', main(2).train_err, main(2).test_err);

fprintf('sub:\ttrain dist: %f\ttest dist: %f\n', sub(1).train_err, sub(1).test_err);
fprintf('sub:\ttrain err: %f\ttest err: %f\n', sub(2).train_err, sub(2).test_err);

fprintf('Xsub:\ttrain dist: %f\ttest dist: %f\n', sub_only(1).train_err, sub_only(1).test_err);
fprintf('Xsub:\ttrain err: %f\ttest err: %f\n', sub_only(2).train_err, sub_only(2).test_err);

%%
function [main, sub] =...
    evaluate_decoder(pre, train, test, X, ks, subset, err_funcs, func_names, train_frac)
train_slice = randperm(length(ks))/length(ks) < train_frac;
X = pre(X);
X_train = X(train_slice,:); X_test = X(~train_slice,:);
ks_train = ks(train_slice); ks_test = ks(~train_slice);
model = train(X_train, ks_train);
pred_train = test(model, X_train);
pred_test = test(model, X_test);

for i = 1:length(err_funcs)
    err = err_funcs{i};
    main(i).err_type = func_names{i};
    main(i).train_err = err(ks_train, pred_train);
    main(i).test_err = err(ks_test, pred_test);

    if ~isempty(subset)
        subset_train = subset(train_slice);
        subset_test = subset(~train_slice);
        sub(i).err_type = func_names{i};
        sub(i).train_err = err(ks_train(subset_train), pred_train(subset_train));
        sub(i).test_err = err(ks_test(subset_test), pred_test(subset_test));
    end
end
end