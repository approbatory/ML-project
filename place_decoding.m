clear;

shuf_train = @(tr) @(X,k) tr(shuffle(X,k),k);
labelshuf_train = @(tr) @(X,k) tr(X, k(randperm(length(k))));
shuf_alg = @(a) struct('name', [a.name ' shuf'], 'pre', a.pre, 'train', shuf_train(a.train), 'test', a.test);
labelshuf_alg = @(a) struct('name', [a.name ' labelshuf'], 'pre', a.pre, 'train', labelshuf_train(a.train), 'test', a.test);

mv_pre = @(X) int64(X ~= 0);
mv_train = @(X_train, ks_train) fitcnb(X_train, ks_train, 'Distribution', 'mvmn',...
    'CategoricalPredictors', 'all');
mv_test = @(model, X_test) predict(model, X_test);

ecoc_pre = @(X) X;
ecoc_pre_binarized = @(X) double(X ~= 0);
t = templateSVM('Standardize',1,'KernelFunction','linear');
t2 = templateSVM('Standardize',1,'KernelFunction','polynomial', 'PolynomialOrder', 2);
ecoc_train = @(X_train, ks_train) fitcecoc(X_train, ks_train, 'Learners', t);
ecoc_train2 = @(X_train, ks_train) fitcecoc(X_train, ks_train, 'Learners', t2);
ecoc_test = @(model, X_test) predict(model, X_test);

mvnb = struct('name', 'Bernoulli NB', 'pre', mv_pre, 'train', mv_train, 'test', mv_test);
ecoc = struct('name', 'ECOC SVM', 'pre', ecoc_pre, 'train', ecoc_train, 'test', ecoc_test);
ecoc2 = struct('name', 'ECOC SVM quadratic', 'pre', ecoc_pre, 'train', ecoc_train2, 'test', ecoc_test);
ecoc_binarized = struct('name', 'ECOC SVM binarized\_input', 'pre', ecoc_pre_binarized, 'train', ecoc_train, 'test', ecoc_test);

algs = [mvnb, shuf_alg(mvnb), ecoc, shuf_alg(ecoc)];%, ecoc_binarized, shuf_alg(ecoc_binarized)];%, ecoc2, shuf_alg(ecoc2)];
%%
directory = '../c14m4';

labels{1} = 'd15, ego left';
days{1} = 'c14m4d15';

labels{2} = 'd16, ego left to allo south';
days{2} = 'c14m4d16';

labels{3} = 'd17, allo south';
days{3} = 'c14m4d17';

par_loops = 16;

for ix = 1:length(days)
    disp(labels{ix});
    decoding_results(ix).label = labels{ix};
    ds = quick_ds(fullfile(directory,days{ix}), 'deprobe', 'nocells');
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
    
    
    
    X = cell2mat(gen_place_decoding_X(ds));
    [ks, D] = bin_space(cell2mat(preprocess_xy(ds)));
    K = size(D,1);
    training_frac = 0.7;
    
    dist_func = @(k,p) D(sub2ind(size(D), k, p));
    err_func = @(k,p) k~=p;
    mean_dist = @(k,p) mean(dist_func(k,p));
    mean_err = @(k,p) mean(err_func(k,p));
    keep_it = @(k,p) [k,p];
    
    err_funcs = {mean_dist, mean_err};%, @confusionmat, keep_it};
    func_names = {'dist', '0-1'};%, 'conf', 'keep'};
    
    for i = 1:length(algs)
        fprintf('Running %s\n', algs(i).name);
        pre = algs(i).pre; train = algs(i).train; test = algs(i).test;
        [main, sub] = ...
            evaluate_decoder(pre, train, test, X, ks, eval_over,...
            err_funcs, func_names, training_frac, par_loops);
        sub_only = evaluate_decoder(pre, train, test,...
            X(eval_over,:), ks(eval_over), [],...
            err_funcs, func_names, training_frac, par_loops);
        decoding_results(ix).res(i).alg = algs(i);
        decoding_results(ix).res(i).division(1).desc = 'full';
        decoding_results(ix).res(i).division(1).out = main;
        decoding_results(ix).res(i).division(2).desc = 'full, test moving';
        decoding_results(ix).res(i).division(2).out = sub;
        decoding_results(ix).res(i).division(3).desc = 'moving';
        decoding_results(ix).res(i).division(3).out = sub_only;
        %fprintf('main:\ttrain dist: %f\ttest dist: %f\n', main(1).train_err, main(1).test_err);
        %fprintf('main:\ttrain err: %f\ttest err: %f\n', main(2).train_err, main(2).test_err);
        
        %fprintf('sub:\ttrain dist: %f\ttest dist: %f\n', sub(1).train_err, sub(1).test_err);
        %fprintf('sub:\ttrain err: %f\ttest err: %f\n', sub(2).train_err, sub(2).test_err);
        
        %fprintf('Xsub:\ttrain dist: %f\ttest dist: %f\n', sub_only(1).train_err, sub_only(1).test_err);
        %fprintf('Xsub:\ttrain err: %f\ttest err: %f\n', sub_only(2).train_err, sub_only(2).test_err);
        fprintf('\n');
    end
end
save multiple_decoding_results.mat decoding_results;
%%
function [main, sub] =...
    evaluate_decoder(pre, train, test, X, ks, subset, err_funcs, func_names, train_frac, par_loops)
if ~exist('par_loops', 'var')
    par_loops = 1;
end
X = pre(X);
tic
parfor par_ix = 1:par_loops
    train_slice{par_ix} = randperm(length(ks))/length(ks) < train_frac;
    X_train = X(train_slice{par_ix},:); X_test = X(~train_slice{par_ix},:);
    ks_train{par_ix} = ks(train_slice{par_ix}); ks_test{par_ix} = ks(~train_slice{par_ix});
    model = train(X_train, ks_train{par_ix});
    pred_train{par_ix} = test(model, X_train);
    pred_test{par_ix} = test(model, X_test);
end
toc
for par_ix = 1:par_loops
    for i = 1:length(err_funcs)
        err = err_funcs{i};
        main(i).err_type = func_names{i};
        main(i).train_err{par_ix} = err(ks_train{par_ix}, pred_train{par_ix});
        main(i).test_err{par_ix} = err(ks_test{par_ix}, pred_test{par_ix});
        
        if ~isempty(subset)
            subset_train = subset(train_slice{par_ix});
            subset_test = subset(~train_slice{par_ix});
            sub(i).err_type = func_names{i};
            sub(i).train_err{par_ix} = err(ks_train{par_ix}(subset_train), pred_train{par_ix}(subset_train));
            sub(i).test_err{par_ix} = err(ks_test{par_ix}(subset_test), pred_test{par_ix}(subset_test));
        end
    end
end
end
