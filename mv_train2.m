function model = mv_train2(X_train, ks_train)
ks_train = ks_train(:)';
K_vals = unique(ks_train);
model.log_prior = sum(ks_train' == K_vals) + 1;
model.log_prior = log(model.log_prior / sum(model.log_prior));
model.log_conditional = zeros(length(K_vals), size(X_train,2));
for i = 1:length(K_vals)
    X = X_train(ks_train == K_vals(i), :);
    model.log_conditional(i,:) = (sum(X) + 1)./(size(X,1) + 2);
end
model.log_flip_conditional = log(1 - model.log_conditional);
model.log_conditional = log(model.log_conditional);
model.K_vals = K_vals;
end