%% declarations
clear;
SHOW_ERRMAP = false;
SHOW_ERR = true;
Lambda = 0.1; %lambda value chosen by cross-validation
Lambda_svm = 0.2; %chosen for SVM
%pre = @(X) double(X ~= 0);
pre = @(X) X;
%train = @(X_train, ks_train) fitclinear(X_train, ks_train,...
%    'Learner', 'logistic', 'ClassNames', [1 2],...
%    'Regularization', 'lasso', 'Lambda', Lambda,...
%    'IterationLimit', 1000, 'Solver', 'sparsa');

%train = @(X_train, ks_train) fitclinear(X_train, ks_train,...
%    'Learner', 'svm', 'ClassNames', [1 2],...
%    'Regularization', 'lasso', 'Lambda', Lambda_svm,...
%    'IterationLimit', 1000, 'Solver', 'sparsa');

train = @(X_train, ks_train) fitcsvm(X_train, ks_train,...
    'KernelFunction', 'polynomial',...
    'PolynomialOrder', 2, 'Standardize', true,...
    'IterationLimit', 10000, 'BoxConstraint', 1,...
    'ClassNames', [1 2], 'Verbose', 0); NONLIN = true;
test = @(model, X_test) predict(model, X_test);

poss = linspace(0,0.4,20); %maybe change

%alg_label = 'logisticL1';
%alg_label = 'svmL1';
alg_label = 'svmP2';
num_runs = 3;
day_label = 'EL2AS';
specific_arm = ', changing arm';

labeling{1} = alg_label;
labeling{2} = [alg_label ' shuf'];
labeling{3} = [alg_label ' labelshuf'];
specific_arm_labeling = cell(1,num_runs);
for i = 1:num_runs
    specific_arm_labeling{i} = [labeling{i} specific_arm];
end

%% Run training alg
ds = quick_ds(fullfile('../c14m4', 'c14m4d16'), 'deprobe', 'nocells');
to_shuf = [0 1 0];
to_labelshuf = [0 0 1];

ks_end = classify_labels({ds.trials.end});
ks_start = classify_labels({ds.trials.start});
X_all = pre(gen_all_X_at_pos_closest(ds, poss));
[X_part, ks_part, vals, masks] = partition_data(X_all, ks_end, ks_start);
err = cell(num_runs, length(vals));
err_map = cell(num_runs, length(vals));
models = cell(num_runs, length(vals));
%fitinfos = cell(num_runs, length(vals));
for j = 1:length(vals)
    if j == 1
        continue;
    end
    X_all = X_part{j};
    ks = ks_part{j};
    for i = 1:num_runs
        err{i,j} = ones(1,length(poss));
        err_map{i,j} = ones(length(ks),length(poss));
        fprintf('\nj%d i%d: %d %d\n',length(vals), num_runs, j,i);
        parfor k = 1:length(poss)
%        for k = 1:length(poss)
            X = X_all(:,:,k);
            %[my_err(k), my_err_map(:,k), my_models{k}, my_fitinfos{k}]...
            [my_err(k), my_err_map(:,k), my_models{k}]...
                = leave_1_out(X, ks, train, test, to_shuf(i), to_labelshuf(i));
            fprintf('*');
        end
        err{i,j} = my_err; err_map{i,j} = my_err_map;
        models{i,j} = my_models; %fitinfos{i,j} = my_fitinfos;
    end
end

%parfor i = 1:num_runs
%    [poss{i}, err{i}, err_map{i}, vals{i}, masks{i}] = decode_end_alg(...
%        pre, train, test, ds, 0.08, 0.4, to_shuf(i), to_labelshuf(i));
%end
%% Plot error maps
if ~exist('SHOW_ERRMAP', 'var') || SHOW_ERRMAP
    for i = 1:num_runs
        view_errmap(ds, poss, err_map(i,2), vals, masks,...
            day_label, labeling{i});
    end
end
%% plot error as function of position
if ~exist('SHOW_ERR', 'var') || SHOW_ERR
    figure;
    ylim([0 inf]);
    hold on;
    grid on;
    formatting = {'-x' '-or' '-og'};
    for i = 1:num_runs
        plot(poss, err{i,2}, formatting{i}, 'DisplayName', specific_arm_labeling{i});
    end
    
    xlabel('arm position');
    ylabel('err');
    title([alg_label ' errs']);
    legend(gca, 'show', 'Location', 'best');
end
%%
if ~exist('NONLIN', 'var') || ~NONLIN
CV_set = models{1,2};
mean_beta = cell(1,length(CV_set));
for k = 1:length(CV_set)
    mean_beta{k} = 0;
    for i = 1:length(CV_set{k})
        mean_beta{k} = mean_beta{k} + CV_set{k}{i}.Beta;
    end
    mean_beta{k} = mean_beta{k}/length(CV_set{k});
end
mean_beta = cell2mat(mean_beta);
[~,order] = sort(mean_beta(:,end));
mean_beta = mean_beta(order,:);
figure; imagesc(mean_beta);
%% Looking at fit params
CV_set = models{1,2}{end};
final_err = err{1,2}(end);
betas = zeros(ds.num_cells, length(CV_set));
for i = 1:length(CV_set)
    betas(:,i) = CV_set{i}.Beta;
end
used_features = find(any(betas~=0,2))';
fprintf('Lambda_svm: %f\t', Lambda_svm);
fprintf('nonzero betas: %d\t', sum(any(betas~=0,2)));
fprintf('residual error: %f\t', final_err);
fprintf('used features: '); disp(used_features);

%looking at correlations:
r_neurons = corr(ks_part{2}', X_part{2}(:,:,end));
[~, order] = sort(r_neurons.^2, 'descend');
clean_order = order(~isnan(r_neurons(order)));

[~, feature_order] = sort(r_neurons(used_features).^2, 'descend');
ranked_features = used_features(feature_order);
end