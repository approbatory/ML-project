function model = pnb_train(X, ks)

[~, n_pred] = size(X);

[K_set, ~, k_ints] = unique(ks);
n_cat = numel(K_set);

mean_resp = zeros(n_cat, n_pred);
cat_frac = zeros(n_cat, 1);

for k_i = 1:n_cat
    k_filt = k_ints == k_i;
    Xk = X(k_filt, :);
    mean_resp(k_i,:) = mean(Xk,1) + 1e-9;
    cat_frac(k_i) = mean(k_filt);
end

model.mean_resp = mean_resp;
model.cat_frac = cat_frac;
model.K_set = K_set;
