function dprime_multi(X, ks)

K = unique(ks);
n_categories = numel(K);
X_class = cell(1,n_categories);
for k_i = 1:n_categories
    category = K(k_i);
    X_class{k_i} = X(ks == category,:);
end

end




function dp2 = dprime(X_a, X_b)
[M_a, N_a] = size(X_a);
[M_b, N_b] = size(X_b);
assert(N_a==N_b, 'Must have same number of predictors %d ~= %d', N_a, N_b);
N = N_a;

ks = [zeros(1,M_a) ones(1,M_b)];
X = [X_a ; X_b];
alg = my_algs('lda');
model = alg.train(X, ks);
dp2 = model.mahal(model.Mu);


if false
    mu_a = mean(X_a);
    mu_b = mean(X_b);
    delta_mu = mu_a - mu_b;
    Sigma2_a = cov(X_a);
    Sigma2_b = cov(X_b);
    Sigma2 = (Sigma2_a + Sigma2_b)/2;
    w_opt = Sigma2\delta_mu;
    w_diag = diag(Sigma2)\delta_mu;
    dp2 = delta_mu.' * w_opt;
    dp_shuffle2 = delta_mu.' * w_diag;
    dp_diagonal2 = dp2^2 / (w_diag.' * Sigma2 * w_diag);
end
end