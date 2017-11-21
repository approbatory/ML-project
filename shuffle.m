function X_out = shuffle(X_in, ks)
    X_out = X_in;
    unique_ks = unique(ks);
    for cell_num = 1:size(X_in, 1)
        for k = unique_ks
            indexes = find(ks == k);
            X_out(cell_num, indexes) = X_in(cell_num, indexes(randperm(size(indexes, 2))));
        end
    end
end



%%Ideas for improving shuffle:
% For the cross-validation, the training set must be shuffled each time,
% but this shuffle can be reduced by only reshuffling the class that the
% excluded training example belongs to. Perhaps shuffling the entire set is
% better for getting a statistical average