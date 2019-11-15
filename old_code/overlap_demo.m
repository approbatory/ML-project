M = 2000;
scale = 10;
S = 7;

X1 = [randn(M,1), randn(M,1)*scale];
X1 = X1 * [1 -1; 1 1]/sqrt(2) + [S/2, S/2];

X2 = [randn(M,1), randn(M,1)*scale];
X2 = X2 * [1 -1; 1 1]/sqrt(2) - [S/2, S/2];

X = [X1; X2];
ks1 = zeros(M,1); ks2 = ones(M,1);
ks = [ks1; ks2];

p = randperm(size(X,1));
X = X(p,:);
ks = ks(p);
svm = my_algs('linsvm');
svm_shuf = my_algs('linsvm', 'shuf');
[tr, te] = evaluate_alg(svm, X, ks, 'eval_f', @(k,p) mean(k(:)~=p(:)), 'par_loops', M);
model = svm.train(X,ks);
[tr_shuf, te_shuf] = evaluate_alg(svm_shuf, X, ks, 'eval_f', @(k,p) mean(k(:)~=p(:)), 'par_loops', M);
model_shuf = svm_shuf.train(X,ks);

figure
subplot(2,1,1);
plot(X1(:,1), X1(:,2), '.b');
hold on
plot(X2(:,1), X2(:,2), '.r');
axis equal
xl = xlim; yl = ylim;
if model.Beta(2) == 0
    h = line([-model.Bias/model.Beta(1), -model.Bias/model.Beta(1)], ylim);
else
    h = refline(get_m(model),get_b(model));
end
h.Color = 'g';
h.LineWidth = 2;
title(sprintf('original, test error: %.2f \\pm %.2f %%, shuffled model: %.2f \\pm %.2f %%', mean(te)*100, std(te)*100, mean(te_shuf)*100, std(te_shuf)*100));
xlim(xl); ylim(yl);

xl = xlim; yl = ylim;
if model_shuf.Beta(2) == 0
    h = line([-model_shuf.Bias/model_shuf.Beta(1), -model_shuf.Bias/model_shuf.Beta(1)], ylim);
else
    h = refline(get_m(model_shuf),get_b(model_shuf));
end
h.Color = 'k';
h.LineWidth = 1;
xlim(xl); ylim(yl);

legend 0 1 location best

X1_shuf = shuffle(X1, ks1);
X2_shuf = shuffle(X2, ks2);
X_shuf = shuffle(X, ks);

[tr_all_shuf, te_all_shuf] = evaluate_alg(svm, X_shuf, ks, 'eval_f', @(k,p) mean(k(:)~=p(:)), 'par_loops', M);

subplot(2,1,2);
plot(X1_shuf(:,1), X1_shuf(:,2), '.b');
hold on
plot(X2_shuf(:,1), X2_shuf(:,2), '.r');
axis equal
xl = xlim; yl = ylim;
if model_shuf.Beta(2) == 0
    h = line([-model_shuf.Bias/model_shuf.Beta(1), -model_shuf.Bias/model_shuf.Beta(1)], ylim);
else
    h = refline(get_m(model_shuf),get_b(model_shuf));
end
h.Color = 'k';
h.LineWidth = 2;
legend 0 1 location best
title(sprintf('intra-class shuffled, test error: %.2f \\pm %.2f %%', mean(te_all_shuf)*100, std(te_all_shuf)*100));
xlim(xl); ylim(yl);


function m = get_m(model)
d = model.Beta(2);
if model.Beta(2) == 0
    d = eps;
end
m = -model.Beta(1)/d;
end

function b = get_b(model)
d = model.Beta(2);
if model.Beta(2) == 0
    d = eps;
end
b = -model.Bias/d;
end