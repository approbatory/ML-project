cohort11 = auto_dayset({'c11m1', 'c11m2', 'c11m3', 'c11m5'});
%c11m1d13

day = cohort11{1}(2);
direc = pwd;
cd(fullfile(day.directory, day.day));
ds = DaySummary(data_sources, 'cm01-fix');
cd(direc);

%%
traces = cell2mat({ds.trials.traces});
traces_at_gate_open = traces(:, ds.trial_indices(:,2));

%trial mask
mask = (strcmp({ds.trials.start}, 'west') & strcmp({ds.trials.turn}, 'right'));
X = traces_at_gate_open(:, mask)';

r = [ds.trials.correct];
ks = r(mask);

%%
reps = 2048;
train_frac = 0.7;
[M,N] = size(X);
parfor i = 1:reps
    partition = randperm(M) <= (train_frac*M);
    X_train = X(partition,:); X_test = X(~partition,:);
    ks_train = ks(partition); ks_test = ks(~partition);
    ks_scrambled = ks(randperm(length(ks)));
    ks_scrambled_train = ks_scrambled(partition); ks_scrambled_test = ks_scrambled(~partition);
    model = fitclinear(X_train, ks_train, 'Learner', 'svm', 'Solver', 'sparsa');
    model_scrambled = fitclinear(X_train, ks_scrambled_train, 'Learner', 'svm', 'Solver', 'sparsa');
    ps_train = predict(model, X_train);
    ps_test = predict(model, X_test);
    train_err(i) = mean(ps_train(:)~=ks_train(:));
    test_err(i) = mean(ps_test(:)~=ks_test(:));
    ps_scrambled_train = predict(model_scrambled, X_train);
    ps_scrambed_test = predict(model_scrambled, X_test);
    train_scrambled_err(i) = mean(ps_scrambled_train(:)~=ks_scrambled_train(:));
    test_scrambled_err(i) = mean(ps_scrambed_test(:)~=ks_scrambled_test(:));
end
fprintf('Decoding correctness (reward/error), at t = gate open, using fluorescence traces:\n');
fprintf('\t Only from trials starting at west (changing arm) and turning right (going south)\n');
fprintf('Train error: %f +- %f (%f +- %f scrambled)\n', mean(train_err),...
    std(train_err)/sqrt(length(train_err)), mean(train_scrambled_err), std(train_scrambled_err)/sqrt(length(train_scrambled_err)));
fprintf(' Test error: %f +- %f (%f +- %f scrambled)\n', mean(test_err), std(test_err)/sqrt(length(test_err)),...
    mean(test_scrambled_err), std(test_scrambled_err)/sqrt(length(test_scrambled_err)));
