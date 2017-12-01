%function [err, err_map, models, fitinfos] = leave_1_out(X, ks, train, test,...
%    do_shuffle, do_label_shuffle)
function [err, err_map, models] = leave_1_out(X, ks, train, test,...
    do_shuffle, do_label_shuffle)
M = length(ks);
err = 0;
err_map = zeros(M,1);
models = cell(1,M);
%fitinfos = cell(1,M);
mask = false(1,M);
count = 0;
for i = 1:M
    mask(i) = true;
    
    X_train = X(~mask,:);
    ks_train = ks(~mask);
    if do_shuffle
        X_train = shuffle(X_train, ks_train);
    end
    if do_label_shuffle
        ks_train = ks_train(randperm(length(ks_train)));
    end
    
    X_test = X(mask,:);
    ks_test = ks(mask);
    
%    [models{i}, fitinfos{i}] = train(X_train, ks_train);
    models{i} = train(X_train, ks_train);
    ks_predicted = test(models{i}, X_test);
    mistakes = sum(ks_test ~= ks_predicted);
    
    total = length(ks_test);
    count = count + 1;
    err = err + mistakes/total;
    err_map(i) = mistakes;
    
    mask(i) = false;
end
err = err/count;
end