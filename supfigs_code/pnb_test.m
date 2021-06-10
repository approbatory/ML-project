function [ps, ll] = pnb_test(model, X)
[n_obs, ~] = size(X);
n_cat = numel(model.K_set);

ll = zeros(n_obs, n_cat);

for k_i = 1:n_cat
    lprior = log(model.cat_frac(k_i));
    q = sum(X .* log(model.mean_resp(k_i,:)), 2) - sum(model.mean_resp(k_i,:));
    ll(:,k_i) = q + lprior;
end

[~, ps] = max(ll, [], 2);
ps = model.K_set(ps);
