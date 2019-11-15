function [one_v_all, one_v_all_d, one_v_one, one_v_one_d] = dprime_multi(X, ks)

K = unique(ks);
n_categories = numel(K);
X_class = cell(n_categories,1);
for k_i = 1:n_categories
    category = K(k_i);
    X_class{k_i} = X(ks == category,:);
end

one_v_all = zeros(1,n_categories);
one_v_all_d = zeros(1,n_categories);
for k_i = 1:n_categories
    [one_v_all(k_i), one_v_all_d(k_i)] = dprime(X_class{k_i}, cell2mat(X_class((1:n_categories) ~= k_i)));
end

one_v_one = zeros(n_categories);
one_v_one_d = zeros(n_categories);
for k_i = 1:n_categories
    for k_j = 1:n_categories
        if k_i == k_j
             one_v_one(k_i, k_j) = Inf;
             one_v_one_d(k_i, k_j) = Inf;
        else
            [one_v_one(k_i, k_j), one_v_one_d(k_i, k_j)] = dprime(X_class{k_i}, X_class{k_j});
        end
    end
end
one_v_one_closest = min(one_v_one);
end




function [dp2, dps2] = dprime(X_a, X_b)
[M_a, N_a] = size(X_a);
[M_b, N_b] = size(X_b);
assert(N_a==N_b, 'Must have same number of predictors %d ~= %d', N_a, N_b);

if false
    ks = [zeros(1,M_a) ones(1,M_b)];
    X = [X_a ; X_b];
    alg = my_algs('lda');
    model = alg.train(X, ks);
    dp2 = model.mahal(model.Mu);
    dp2 = dp2(2);
    
    alg = my_algs('dlda');
    model = alg.train(X, ks);
    dps2 = model.mahal(model.Mu);
    dps2 = dps2(2);
else
    mu_a = mean(X_a).';
    mu_b = mean(X_b).';
    delta_mu = mu_a - mu_b;
    %Sigma2_a = cov(X_a);
    %Sigma2_b = cov(X_b);
    %Sigma2 = (Sigma2_a + Sigma2_b)/2;
    Sigma2 = cov([(X_a - mu_a.'); (X_b - mu_b.')]);
    w_opt = Sigma2\delta_mu;
    w_diag = delta_mu./diag(Sigma2);
    dp2 = delta_mu.' * w_opt;
    dp_shuffle2 = delta_mu.' * w_diag;
    dp_diagonal2 = dp2^2 / (w_diag.' * Sigma2 * w_diag);
    
    dps2 = dp_shuffle2;
end
end