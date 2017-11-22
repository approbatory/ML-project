%running Multinomial_NB
%run_alg(@mnb_train, @mnb_test, 'MultinomialNB', 'step', 0.04);

%%%%%%running GaussNB <- can't run this because some features have 0 variance
%another Multinomial_NB, arbitrarily discretize it, using a factor of 40
%(approx. value of mn dist. n parameter)
m2nb_pre = @(X) round(40*X);
m2nb_train = @(X_train, ks_train) fitcnb(X_train, ks_train, 'Distribution', 'mn');
m2nb_test = @(model, X_test) predict(model, X_test);
run_alg(m2nb_train, m2nb_test, 'MultinomialNB2', 'preprocessor', @round, 'step', 0.04);

%running SVM
svm_train = @(X_train, ks_train) fitclinear(X_train, ks_train,...
    'Learner', 'svm', 'ClassNames', [1 2]);
svm_test = @(model, X_test) predict(model, X_test);
run_alg(svm_train, svm_test, 'SVM', 'step', 0.04);
disp('NOW SHUFFLING');
run_alg(svm_train, svm_test, 'SVM', 'shuffle', 'step', 0.04);



%-----------------------
function params = mnb_train(X_train, ks_train)
    [params.log_prior, params.log_conditional] =...
        multinom_nb_encode(X_train, ks_train, numel(unique(ks_train)));
end
function pred = mnb_test(params, X_test)
    pred = multinom_nb_decode(X_test, params.log_prior,...
        params.log_conditional);
end