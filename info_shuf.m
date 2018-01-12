directory = '../c14m4';

labels{1} = 'd15, ego left';
days{1} = 'c14m4d15';

labels{2} = 'd16, ego left to allo south';
days{2} = 'c14m4d16';

labels{3} = 'd17, allo south';
days{3} = 'c14m4d17';


tlin = templateLinear('Learner', 'svm',...
    'Regularization', 'lasso', 'Lambda', 0.002,...
    'Solver', 'sparsa');
ecoc_trainlin = @(X_train, ks_train) fitcecoc(X_train, ks_train, 'Learners', tlin);
ecoclin = struct('name', 'ECOC SVM lin',...
                 'pre', @(X) X,...
                 'train', ecoc_trainlin,...
                 'test', @predict);
alg = ecoclin;
train_frac = 0.7;

process_days = 1;
for ix = process_days
    disp(labels{ix});
    ds = quick_ds(fullfile(directory,days{ix}), 'deprobe', 'nocells');
    X = sparse(cell2mat(gen_place_decoding_X(ds)));
    [M,N] = size(X);
    xy = cell2mat({ds.trials.centroids}');
    Xb = X ~= 0;
    [ks, D] = bin_space(cell2mat(preprocess_xy(ds)));
    err_func = @(k,p) mean(D(sub2ind(size(D),k,p)));
    K = size(D,1);
    col_muti = cellfun(@(x) muti(ks,x), num2cell(Xb,1));
    %[~, order] = sort(col_muti);
    %X_ord = X(:,order);
    %neu = full(X_ord(:,end));
    %vis(neu~=0, xy);
    perf = @(a) performance(a, X, ks, err_func, train_frac);
    
    [base_train_err, base_test_err] = perf(alg);
    fprintf('Baseline:\ttr: %f\tte: %f\n', base_train_err, base_test_err);
    
    totshuf_alg = modify_alg(alg, @shuffle, ' totshuf');
    [totshuf_train_err, totshuf_test_err] = perf(totshuf_alg);
    fprintf('Total shuffle:\ttr: %f\tte: %f\n', totshuf_train_err, totshuf_test_err);
    singleshuf_train_err = zeros(1,N); singleshuf_test_err = zeros(1,N);
    parfor j = 1:N
        %featmask = false(1,N);
        %featmask(j) = true;
        featmask = randperm(N) <= j;
        singleshuf_alg = modify_alg(alg, @(X,ks) shufgen(X,ks,featmask), sprintf(' shuf%d',j));
        [singleshuf_train_err(j), singleshuf_test_err(j)] = perf(singleshuf_alg);
        fprintf('Done j=%d of %d\n', j, N);
    end
    change_in_train_err = singleshuf_train_err - base_train_err; 
    change_in_test_err = singleshuf_test_err - base_test_err;
end


function mod_alg = modify_alg(alg, modder, modlabel)
mod_alg = struct('name', [alg.name modlabel],...
                 'pre', alg.pre,...
                 'train', @(X,ks) alg.train(modder(X,ks),ks),...
                 'test', alg.test);
end


function [train_err, test_err] = performance(alg, X, ks, func, train_frac)
X = alg.pre(X);
train_slice = randperm(length(ks))/length(ks) < train_frac;
X_train = X(train_slice,:); X_test = X(~train_slice,:);
ks_train = ks(train_slice); ks_test = ks(~train_slice);
model = alg.train(X_train, ks_train);
pred_train = alg.test(model, X_train);
pred_test = alg.test(model, X_test);

train_err = func(ks_train, pred_train);
test_err = func(ks_test, pred_test);
end


function vis(neu, xy)
xy_on = xy(neu,:);
xy_off = xy(~neu,:);
alpha = 0.1;
figure; scatter(xy_off(:,1), xy_off(:,2), 1,...
    'MarkerEdgeColor', 'k', 'MarkerEdgeAlpha', alpha,...
    'MarkerFaceColor', 'k', 'MarkerFaceAlpha', alpha);
hold on; scatter(xy_on(:,1), xy_on(:,2), 1,...
    'MarkerEdgeColor', 'r', 'MarkerEdgeAlpha', alpha,...
    'MarkerFaceColor', 'r', 'MarkerFaceAlpha', alpha);
hold off;
end