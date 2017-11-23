function model = mv_train2(X_train, ks_train)
K_vals = unique(ks_train);
model.log_prior = sum(ks_train' == K_vals) + 1;
model.log_prior = log(model.log_prior / sum(model.log_prior));
model.log_conditional = zeros(length(K_vals), size(X_train,2));
for k = K_vals
    X = X_train(ks_train == k, :);
    model.log_conditional(k,:) = log((sum(X) + 1)./(size(X,1) + 2));
end
end