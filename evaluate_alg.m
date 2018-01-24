function [model, train_err, test_err] = evaluate_alg(alg, X, ks, eval_f, train_frac, par_loops)
if ~exist('par_loops', 'var')
    par_loops = 1;
end
tic
X = alg.pre(X);
parfor par_ix = 1:par_loops
    train_slice = randperm(length(ks))/length(ks) < train_frac;
    X_train = X(train_slice,:); X_test = X(~train_slice,:);
    ks_train = ks(train_slice); ks_test = ks(~train_slice);
    model{par_ix} = train(X_train, ks_train);
    pred_train = test(model{par_ix}, X_train);
    pred_test = test(model{par_ix}, X_test);
    train_err{par_ix} = eval_f(ks_train, pred_train);
    test_err{par_ix} = eval_f(ks_test, pred_test);
end
toc
end