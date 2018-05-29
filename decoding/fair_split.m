function subset = fair_split(ks, train_frac)
subset = false(size(ks));
K = unique(ks(:))';
for k = K
    num_k = sum(ks == k);
    sel = randperm(num_k)/num_k <= train_frac;
    subset(ks == k) = sel;
end
end