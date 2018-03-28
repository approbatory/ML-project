function model = mn_train2(X_train,ks_train)
[M,N] = size(X_train);
K_vals = unique(ks_train);
K = length(K_vals);
model.log_prior = log(sum(K_vals==ks_train')/M);
model.log_conditional = zeros(N,K);
for k = 1:K
    model.log_conditional(:,k) = sum(X_train(ks_train==K_vals(k),:))+1;
end
model.log_conditional = log(model.log_conditional./sum(model.log_conditional));
end

