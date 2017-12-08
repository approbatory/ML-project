function X = shuffle(X, ks)
ks = ks(:)';
K_vals = unique(ks);
N = size(X,2);
for j = 1:N
    for k = K_vals
        is = find(ks == k);
        X(is, j) = X(is(randperm(length(is))),j);
    end
end
end