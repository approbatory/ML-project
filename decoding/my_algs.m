function algs = my_algs(names, settings, all2all, Lambda)
if (nargin == 2) && isnumeric(settings)
    Lambda = settings;
    settings = 'original';
end

%declarations
shuf_train = @(tr) @(X,k) tr(shuffle(X,k),k);
labelshuf_train = @(tr) @(X,k) tr(X, k(randperm(length(k))));
shuf_alg = @(a) struct('name', [a.name ' (shuffled)'], 'pre', a.pre, 'train', shuf_train(a.train), 'test', a.test);
labelshuf_alg = @(a) struct('name', [a.name ' (labelshuf)'], 'pre', a.pre, 'train', labelshuf_train(a.train), 'test', a.test);

mv_pre = @(X) double(X ~= 0);
mv_train = @(X_train, ks_train) fitcnb(X_train, ks_train, 'Distribution', 'mvmn',...
    'CategoricalPredictors', 'all');
mv_test = @(model, X_test) predict(model, X_test);

gnb_pre = @(X)X;
gnb_train = @(X,ks) fitcnb(X,ks);
gnb_test = @predict;

ecoc_pre = @(X) X;
ecoc_pre_binarized = @(X) double(X ~= 0);
t = templateSVM('Standardize',1,'KernelFunction','linear');
t2 = templateSVM('Standardize',1,'KernelFunction','polynomial', 'PolynomialOrder', 2);
if ~exist('Lambda', 'var')
    L = 0;
else
    L = Lambda;
end
tlin = templateLinear('Learner', 'svm',...
    'Regularization', 'lasso', 'Lambda', L,...
    'Solver', 'sparsa');

%if ~exist('Lambda', 'var')
%    L = 0.02;
%else
%    L = Lambda;
%end
linsvm_train = @(X_train, ks_train) fitclinear(X_train, ks_train, 'Learner', 'svm',...
    'Regularization', 'lasso', 'Lambda', L,...
    'Solver', 'sparsa');
linsvm_test = @(model, X_test) predict(model, X_test);

ecoc_train = @(X_train, ks_train) fitcecoc(X_train, ks_train, 'Learners', t);
ecoc_train2 = @(X_train, ks_train) fitcecoc(X_train, ks_train, 'Learners', t2);
ecoc_trainlin = @(X_train, ks_train) fitcecoc(X_train, ks_train, 'Learners', tlin);
ecoc_trainlin_onevsall = @(X_train, ks_train) fitcecoc(X_train, ks_train, 'Coding', 'onevsall', 'Learners', tlin);
ecoc_test = @(model, X_test) predict(model, X_test);



package_fun = @(name, trainer) struct('name', name, 'pre', @(X)X, 'train', trainer, 'test', @predict);

gnb = struct('name', 'Gaussian NB', 'pre', gnb_pre, 'train', gnb_train, 'test', gnb_test);
lda = package_fun('LDA', @(X,ks) fitcdiscr(X,ks,'DiscrimType', 'linear'));
qda = package_fun('QDA', @(X,ks) fitcdiscr(X,ks,'DiscrimType', 'quadratic'));
mvnb = struct('name', 'Bernoulli NB', 'pre', mv_pre, 'train', mv_train, 'test', mv_test);
mvnb2= struct('name', 'Bernoulli NB','pre', mv_pre, 'train',@mv_train2,'test',@mv_test2);
ecoc = struct('name', 'ECOC SVM', 'pre', ecoc_pre, 'train', ecoc_train, 'test', ecoc_test);
ecoc2 = struct('name', 'ECOC SVM quadratic', 'pre', ecoc_pre, 'train', ecoc_train2, 'test', ecoc_test);
ecoclin = struct('name', 'pairwise ECOC SVM', 'pre', ecoc_pre, 'train', ecoc_trainlin, 'test', ecoc_test);
ecoclin_onevsall = struct('name', '1 vs. others ECOC SVM', 'pre', ecoc_pre, 'train', ecoc_trainlin_onevsall, 'test', ecoc_test);
ecoc_binarized = struct('name', 'ECOC SVM binarized\_input', 'pre', ecoc_pre_binarized, 'train', ecoc_train, 'test', ecoc_test);
ecoclin_binarized = struct('name', 'ECOC SVM lin binarized\_input', 'pre', ecoc_pre_binarized, 'train', ecoc_trainlin, 'test', ecoc_test);

linsvm = struct('name', 'lin SVM', 'pre', @(x)x, 'train', linsvm_train, 'test', linsvm_test);

if ~exist('settings', 'var')
    settings = 'original';
end
if ~exist('all2all', 'var')
    all2all = true;
end

if ~iscell(names)
    names = {names};
end
if ~iscell(settings)
    settings = {settings};
end
if numel(names)~=numel(settings) && ~all2all
    error('names and settings must match in length');
end
%lookups
if all2all
    for i = 1:numel(names)
        for j = 1:numel(settings)
            algs(j,i) = create(names{i}, settings{j});
        end
    end
else
for i = 1:numel(names)
    algs(i) = create(names{i}, settings{i});
end
end

    function alg = create(name, setting)
        switch name
            case 'gnb'
                alg = gnb;
            case 'lda'
                alg = lda;
            case 'qda'
                alg = qda;
            case 'mvnb'
                alg = mvnb;
            case 'mvnb2'
                alg = mvnb2;
            case 'ecoc'
                alg = ecoc;
            case 'ecoc2'
                alg = ecoc2;
            case 'ecoclin'
                alg = ecoclin;
            case 'ecoclin_onevsall'
                alg = ecoclin_onevsall;
            case 'ecoc_binarized'
                alg = ecoc_binarized;
            case 'ecoclin_binarized'
                alg = ecoclin_binarized;
            case 'linsvm'
                alg = linsvm;
            otherwise
                error('requested alg name not available: %s', name);
        end
        switch setting
            case 'original'
            case 'shuf'
                alg = shuf_alg(alg);
            case 'labelshuf'
                alg = labelshuf_alg(alg);
            otherwise
                error('requested setting not available: %s', setting);
        end
    end
end



function pred = mv_test2(model, X_test)
[~, pred] = max(model.log_prior' + sum(model.log_flip_conditional,2) +...
    (model.log_conditional - model.log_flip_conditional)*X_test');
pred = pred';
pred = model.K_vals(pred);
end


function model = mv_train2(X_train, ks_train)
ks_train = ks_train(:)';
K_vals = unique(ks_train);
model.log_prior = sum(ks_train' == K_vals) + 1;
model.log_prior = log(model.log_prior / sum(model.log_prior));
model.log_conditional = zeros(length(K_vals), size(X_train,2));
for i = 1:length(K_vals)
    X = X_train(ks_train == K_vals(i), :);
    model.log_conditional(i,:) = (sum(X) + 1)./(size(X,1) + 2);
end
model.log_flip_conditional = log(1 - model.log_conditional);
model.log_conditional = log(model.log_conditional);
model.K_vals = K_vals;
end