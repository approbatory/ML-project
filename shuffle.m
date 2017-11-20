function X_out = shuffle(X_in, ks)
    X_out = X_in;
    for cell_num = 1:size(X_in, 1)
        for k = unique(ks)
            indexes = find(ks == k);
            X_out(cell_num, indexes) = X_in(cell_num, indexes(randperm(size(indexes, 2))));
        end
    end
end
