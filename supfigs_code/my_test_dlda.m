function ps = my_test_dlda(model, X)

X = X ./ model.pred_std;

%norm_mean = permute(model.norm_mean, [3 2 1]);

%[~, ps] = min(sum((X-norm_mean).^2,2),[],3);
[n_obs, ~] = size(X);
dist2 = zeros(n_obs, numel(model.K_set));

for k_i = 1:numel(model.K_set)
    dist2(:, k_i) = sum((X - model.norm_mean(k_i,:)).^2,2);
end
[~, ps] = min(dist2, [], 2);


ps = model.K_set(ps);