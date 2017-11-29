clear;
pre = @(X) X;
%svm_train = @(X_train, ks_train) fitclinear(X_train, ks_train,...
%    'Learner', 'svm', 'ClassNames', [1 2]);
train = @(X_train, ks_train) fitcsvm(X_train, ks_train,...
    'KernelFunction', 'linear', 'Standardize', true,...
    'ClassNames', [1 2]);
test = @(model, X_test) predict(model, X_test);


ds = quick_ds(fullfile('../c14m4', 'c14m4d16'), 'deprobe', 'nocells');
[poss{1}, err{1}, err_map{1}, vals{1}, masks{1}] = decode_end_alg(...
    pre, train, test, ds, 0.08, 0.4, 0);
[poss{2}, err{2}, err_map{2}, vals{2}, masks{2}] = decode_end_alg(...
    pre, train, test, ds, 0.08, 0.4, 1);
[poss{3}, err{3}, err_map{3}, vals{3}, masks{3}] = decode_end_alg(...
    pre, train, test, ds, 0.08, 0.4, 0, 1);

view_errmap(ds, poss{1}, err_map{1}, vals{1}, masks{1},...
    'ego left to allo south', 'SVM');
view_errmap(ds, poss{2}, err_map{2}, vals{2}, masks{2},...
    'ego left to allo south', 'SVM shuf');
view_errmap(ds, poss{3}, err_map{3}, vals{3}, masks{3},...
    'ego left to allo south', 'SVM labelshuf');

figure;
plot(poss{1}, err{1}{2}, '-x', 'DisplayName', 'svm, changing arm');
ylim([0 inf]);
hold on;
plot(poss{2}, err{2}{2}, '-or', 'DisplayName', 'svm shuf, changing arm');
plot(poss{3}, err{3}{2}, '-og', 'DisplayName', 'svm labelshuf, changing arm');

xlabel('arm position');
ylabel('err');
title('SVM errs');
legend(gca, 'show');