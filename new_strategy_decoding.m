%[dayset, matching] = my_daysets('c11m1'); end_arm = 'north'; %mPFC
[dayset, matching] = my_daysets('c14m4'); end_arm = 'south'; %HPC

points = 0:0.1:1;
for i = 1:numel(dayset)
    ds(i) = quick_ds(fullfile(dayset(i).directory, dayset(i).day), 'deprobe', 'nocells');
    if i == 1
        t = zeros(1,ds(i).num_trials);
    elseif i == 3
        t = ones(1,ds(i).num_trials);
    elseif i == 2
        t = [zeros(1,50) ones(1,ds(i).num_trials-50)];
    end
    for j = 1:length(points)
        [X{i,j}, ks{i,j}, ~] = ds_dataset(ds(i), 'filling', 'binary', 'selection', points(j),...
            'trials', strcmp({ds(i).trials.start}, 'east') & strcmp({ds(i).trials.end}, end_arm),...
            'target', num2cell(t));
        X{i,j} = X{i,j}(:, matching(:,i));
    end
end

%% TODO run the decoding with the data
LAMBDA = 0.05; %for HPC

for j = 1:length(points)
    X_train{j} = [X{1,j}; X{3,j}];
    ks_train{j} = cell2mat([ks{1,j}, ks{3,j}]);

    X_test{j} = X{2,j};
    ks_test{j} = cell2mat(ks{2,j});
end

find_errs = @(k,p) k(:) ~= p(:);
%algs = my_algs({'mvnb2', 'linsvm'}, {'original', 'shuf'});
algs(1) = my_algs({'mvnb2'}, {'original'});
algs(1).name = 'Naive Bayes';

svm_train = @(X_train,ks_train) fitclinear(X_train, ks_train, 'Learner', 'svm',...
    'Regularization', 'lasso', 'Lambda', LAMBDA,...
    'Solver', 'sparsa');
algs(2) = struct('name', 'SVM', 'train', svm_train, 'test', @predict, 'pre', @(x)x);
algs(3) = struct('name', 'SVM (tr. shuf)', 'train', @(X,k) svm_train(shuffle(X,k),k), 'test', @predict, 'pre', @(x)x);
figure;
for i = 1:numel(algs)
    disp(algs(i).name);
    alg = algs(i);
    for j = 1:length(points)
        disp(points(j));
        model = alg.train(X_train{j}, ks_train{j});
        train_pred{i}(:,j) = alg.test(model, X_train{j});
        train_errors{i}(:,j) = find_errs(alg.test(model, X_train{j}), ks_train{j});
        test_pred{i}(:,j) = alg.test(model, X_test{j});
        test_errors{i}(:,j) = find_errs(alg.test(model, X_test{j}), ks_test{j});
    end
    subplot(2, numel(algs), i);
    imagesc(train_pred{i});
    set(gca, 'XTickLabel', {'.1', '.3', '.5', '.7', '.9'});
    h = refline([0 48.5]); h.Color = 'r'; h.LineWidth = 2;
    title(['Train set, ' alg.name]);
    colormap(gray);
    subplot(2, numel(algs), i+numel(algs));
    imagesc(test_pred{i});
    set(gca, 'XTickLabel', {'.1', '.3', '.5', '.7', '.9'});
    h = refline([0 24.5]); h.Color = 'r'; h.LineWidth = 2;
    title(['Test set, ' alg.name]);
    colormap(gray);
end