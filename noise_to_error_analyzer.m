dps = 0:0.01:10;
res = ones(size(dps));

for i = 1:numel(dps) 
N = 100000;
X1 = randn(N,1);
X2 = randn(N,1) + dps(i);

%figure; hold on;
%histogram(X1);
%histogram(X2);

X = [X1;X2];
y = [zeros(size(X1)); ones(size(X2))];
[X_tr, X_te, y_tr, y_te] = holdout_selector(0.3, X, y);
alg = my_algs('lda');
model = alg.train(X_tr, y_tr);
err_tr = mean(abs(model.predict(X_tr) - y_tr));
err_te = mean(abs(model.predict(X_te) - y_te));
res(i) = err_te;
if mod(i,10) == 0
    fprintf('%d / %d\n', i, numel(dps));
end
end

%%
figure;
semilogy(dps, res);
title 'error by d'''
xlabel 'd'''
ylabel error

%since the error rate between nearest bins is ~0.3, then the d' between
%them should be ~1