clear;
K = 20;
M = 300000;
N = 200;
%%
dayset = my_daysets('c14m4');
[ds, X, ks, errf] = load_day(dayset(1));
nb = my_algs('mvnb2', 'original');
X = nb.pre(X);
%train_set = randperm(length(ks))/length(ks) < 0.7;
%X_tr = X(train_set,:); X_te = X(~train_set,:);
%ks_tr = ks(train_set); ks_te = ks(~train_set);
model = nb.train(X, ks);
P = exp(model.log_conditional);
pr = exp(model.log_prior);

fake_ks = sum(rand(length(ks),1) > cumsum(pr),2)+1;
fake_X = sparse(double(rand(size(X)) < P(fake_ks,:)));
%%

%P = 0.1*abs(randn(K, N));
%ks = randi(K, M, 1);

%X = sparse(double(rand(M,N) < P(ks,:)));

nb = my_algs('mvnb2', 'original');
ecoc = my_algs('ecoclin', 'original');


[train_err_nb, test_err_nb] = evaluate_alg(nb, fake_X, fake_ks,...
    'eval_f', @mean_dist);

[train_err_ecoc, test_err_ecoc] = evaluate_alg(ecoc, fake_X, fake_ks,...
    'eval_f', @mean_dist);

fprintf('NB\nTrain err was\t%f\nTest err was\t%f\n',...
    train_err_nb, test_err_nb);
fprintf('ECOC\nTrain err was\t%f\nTest err was\t%f\n',...
    train_err_ecoc, test_err_ecoc);


function p = isprob(s, K, h)
p = exp(s'*K*s + h*s);
dh = p*s;
dK = p * (s * s');
end