function model = my_dlda(X, ks)

[~, n_pred] = size(X);

[K_set, ~, k_ints] = unique(ks);
n_cat = numel(K_set);

mean_resp = zeros(n_cat, n_pred);

for k_i = 1:n_cat
    mean_resp(k_i,:) = mean(X(k_ints == k_i, :),1);
end

X_mean_only = mean_resp(k_ints,:);
X_noise = X - X_mean_only;

pred_std = std(X_noise, [], 1);

model.pred_std = pred_std;
model.norm_mean = mean_resp ./ pred_std;
model.K_set = K_set;