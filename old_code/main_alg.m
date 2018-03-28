%% MultivariateNB
mv_pre = @(X) int64(X ~= 0);
mv_train = @(X_train, ks_train) fitcnb(X_train, ks_train, 'Distribution', 'mvmn',...
    'CategoricalPredictors', 'all');
mv_test = @(model, X_test) predict(model, X_test);

run_alg(mv_train, mv_test, 'MultivariateNB', 'preprocessor', mv_pre, 'step', 0.04);
run_alg(mv_train, mv_test, 'MultivariateNB', 'preprocessor', mv_pre, 'step', 0.04, 'shuffle');

run_alg(@mv_train2, @mv_test2, 'MultivariateNB2', 'preprocessor', mv_pre, 'step', 0.04);
run_alg(@mv_train2, @mv_test2, 'MultivariateNB2', 'preprocessor', mv_pre, 'step', 0.04, 'shuffle');
%% MultinomialNB
mn_pre = @(X) round(X);
mn_train = @(X_train, ks_train) fitcnb(X_train, ks_train, 'Distribution', 'mn');
mn_test = @(model, X_test) predict(model, X_test);

run_alg(mn_train, mn_test, 'MultinomialNB', 'preprocessor', mn_pre, 'step', 0.04);
run_alg(mn_train, mn_test, 'MultinomialNB', 'preprocessor', mn_pre, 'step', 0.04, 'shuffle');
%%
run_alg(@mn_train2, @mn_test2, 'MultinomialNB2', 'preprocessor', mn_pre, 'step', 0.04);
run_alg(@mn_train2, @mn_test2, 'MultinomialNB2', 'preprocessor', mn_pre, 'step', 0.04, 'shuffle');
%% Linear SVM + shuffle
%svm_pre = @(X) (X-mean(X))./std(X);
svm_pre = @standardize; 
svm_train = @(X_train, ks_train) fitclinear(X_train, ks_train,...
    'Learner', 'svm', 'ClassNames', [1 2]);
svm_test = @(model, X_test) predict(model, X_test);
%%
run_alg(svm_train, svm_test, 'SVM', 'step', 0.04);
run_alg(svm_train, svm_test, 'SVM', 'shuffle', 'step', 0.04);
%%
run_alg(svm_train, svm_test, 'SVMstd', 'preprocessor', svm_pre, 'step', 0.04);
run_alg(svm_train, svm_test, 'SVMstd', 'preprocessor', svm_pre, 'shuffle', 'step', 0.04);

%%
function X = standardize(X)
S = size(X,3);
X = X - mean(X);
sigma = std(X);
featmask = true(1,size(X,2));
for k = 1:S
    featmask = featmask & (sigma(:,:,k)~=0);
end
X = X(:,featmask,:) ./ sigma(:,featmask,:);
end