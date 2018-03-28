mv_pre = @(X) int64(X ~= 0);
mv_train = @(X_train, ks_train) fitcnb(X_train, ks_train, 'Distribution', 'mvmn',...
    'CategoricalPredictors', 'all');
mv_test = @(model, X_test) predict(model, X_test);

%run_alg(mv_train, mv_test, 'MultivariateNB', 'preprocessor', mv_pre,...
%    'step', 0.08, 'nosave', 'nopar', 'errmaps');


run_alg(@mv_train2, @mv_test2, 'MultivariateNB2', 'preprocessor', mv_pre,...
    'step', 0.08, 'nosave', 'nopar', 'errmaps');

%%
svm_train = @(X_train, ks_train) fitclinear(X_train, ks_train,...
    'Learner', 'svm', 'ClassNames', [1 2]);
svm_test = @(model, X_test) predict(model, X_test);

run_alg(svm_train, svm_test, 'SVM', 'step', 0.08,...
    'nosave', 'nopar', 'errmaps');
run_alg(svm_train, svm_test, 'SVM', 'shuffle', 'step', 0.08,...
    'nosave', 'nopar', 'errmaps');

%% only switch day, changing arm
clear;
svm_pre = @(X) X;
svm_train = @(X_train, ks_train) fitclinear(X_train, ks_train,...
    'Learner', 'svm', 'ClassNames', [1 2]);
svm_test = @(model, X_test) predict(model, X_test);


ds = quick_ds(fullfile('../c14m4', 'c14m4d16'), 'deprobe', 'nocells');
[poss{1}, err{1}, err_map{1}, vals{1}, masks{1}] = decode_end_alg(...
    svm_pre, svm_train, svm_test, ds, 0.08, 0.4, 0);
[poss{2}, err{2}, err_map{2}, vals{2}, masks{2}] = decode_end_alg(...
    svm_pre, svm_train, svm_test, ds, 0.08, 0.4, 1);
view_errmap(ds, poss{1}, err_map{1}(2), vals{1}(2), masks{1}(2),...
    'ego left to allo south', 'SVM');
view_errmap(ds, poss{2}, err_map{2}(2), vals{2}(2), masks{2}(2),...
    'ego left to allo south', 'SVM shuf');

figure;
plot(poss{1}, err{1}{2}, '-x', 'DisplayName', 'svm, changing arm');
ylim([0 0.5]);
hold on;
plot(poss{2}, err{2}{2}, '-or', 'DisplayName', 'svm shuf, changing arm');

xlabel('arm position');
ylabel('err');
title('SVM errs');
legend(gca, 'show');