L = 16; N = 1000;
points = linspace(0,L,1000);
%means = L*rand(1,N)/2;
%means2 = L*rand(1,N)/2 + L/2;
means = L*rand(1,N);
means2 = L*rand(1,N);
stds = ones(1,N);
%base_lambda = 0.1; %this was the case for the two place field case
base_lambda = 0.01;
%Lambdas = @(point) base_lambda .* (exp(-(point - means).^2./(2*stds)) + exp(-(point - means2).^2./(2*stds))); %(position, neuron)
Lambdas = @(point) base_lambda .* (exp(-(point - means).^2./(2*stds)) + exp(-(point - means).^2./(2*stds))); %(position, neuron)

times = 1:100000; speed = 0.01;
trajectory = L*sin(speed*times).^2;
%%
rates = Lambdas(trajectory');

X = double(sparse(poissrnd(rates) > 0));
ks = ceil(trajectory);
ks_LR = ks > 8;
errf = @(k,p) mean(abs(k(:)-p(:)));

%%
algs = my_algs({'mvnb2', 'ecoclin'}, {'original', 'shuf'}, true, 0);
PARLOOPS = 4;
train_err = cell(2,numel(algs)); test_err = cell(2,numel(algs));
for i = 1:numel(algs)
    [train_err{1,i}, test_err{1,i}] = evaluate_alg(algs(i),...
        X, ks, 'eval_f', errf,...
        'train_frac', 0.7, 'par_loops', PARLOOPS);
    [train_err{2,i}, test_err{2,i}] = evaluate_alg(algs(i),...
        X, ks_LR, 'eval_f', errf,...
        'train_frac', 0.7, 'par_loops', PARLOOPS);
end
%train_err(2,:) = train_err(1,:);
%test_err(2,:) = test_err(1,:);
figure;
dayset(1).label = '16 bins';
dayset(2).label = 'L/R';

plotmat('Effect of bin size (1 random place field)', train_err, test_err, algs, dayset, 5);

return;
%%
bin_scheme = [
    1 0 1 0 % 0 0 0 0
    1 0 1 1 % 0 0 0 1
    1 0 0 1 % 0 0 1 1
    1 0 0 0 % 0 0 1 0
    1 1 0 0 % 0 1 1 0
    1 1 0 1 % 0 1 1 1
    1 1 1 1 % 0 1 0 1
    1 1 1 0 % 0 1 0 0
    0 1 1 0 % 1 1 0 0
    0 1 1 1 % 1 1 0 1
    0 1 0 1 % 1 1 1 1
    0 1 0 0 % 1 1 1 0
    0 0 0 0 % 1 0 1 0
    0 0 0 1 % 1 0 1 1
    0 0 1 1 % 1 0 0 1
    0 0 1 0 % 1 0 0 0
    ];
bounds = {8, [4,12], [2,6,10,14], [1,3,5,7,9,11,13,15]};
sharpness = 10;
sg = @(x) 1./(1+exp(-sharpness.*x));
sens = @(x,bounds) sum((mod((1:numel(bounds))',2)*2-1).*sg(x - bounds'),1);
%%
N = 500;
W = 1000;
%W = 20;
sparsity = 0.5;
K = floor(N*sparsity);
pos_words = scode(N,K,W);%randi(2,W,N)-1;
neg_words = scode(N,K,W);%randi(2,W,N)-1;
words = [pos_words; neg_words]; ks_words = [ones(1,W), zeros(1,W)]; words = sparse(words);
self_prob = prod(mean(pos_words).* pos_words + (1-mean(pos_words)) .* (1-pos_words),2);
diff_prob = prod(mean(pos_words).* neg_words + (1-mean(pos_words)) .* (1-neg_words),2);
err_rate = mean(self_prob < diff_prob);
perf(words, ks_words);
%%
N = 500;
%W = 2520;
W = 256;
K = 2;
K2 = 8;
sparsity = 0.01;
F = floor(N*sparsity);
words = sparse(scode(N,F,W));
ks_words = mod(1:W,K);
ks_words2 = mod(1:W,K2);
perf(words, ks_words, ks_words2, true, 'K=2', 'K=8');
%%
M = 10000;%M = 10000;
ks = randi(2,1,M) - 1;
X = zeros(M,N);
for i = 1:M
    X(i,:) = genr(ks(i), pos_words, neg_words, 0, K);
end
%X_noise = rand(size(X));
%X(X_noise > 0.2) = 0; % fires only p of the time
%X = sparse(X);
%X = X .* (rand(size(X)) < 0.1);
errf = @(k,p) mean(k(:)~=p(:));
%%
alpha_size = 2000;
pos_word = scode(N,K,1);
neg_word = scode(N,K,1);
alpha_words = scode(N,floor(N/2),alpha_size);
pos_words = xor(pos_word, alpha_words);
neg_words = xor(neg_word, alpha_words);

self_prob = prod(mean(pos_words).* pos_words + (1-mean(pos_words)) .* (1-pos_words),2);
diff_prob = prod(mean(pos_words).* neg_words + (1-mean(pos_words)) .* (1-neg_words),2);
err_rate = mean(self_prob < diff_prob);
words = [pos_words; neg_words]; ks_words = [ones(1,alpha_size), zeros(1,alpha_size)];
M = 10000;
[X, ks] = data_gen(M, pos_words, neg_words);
%% vis
Y = tsne(full(double(X)), 'Distance', 'hamming');
figure; gscatter(Y(:,1), Y(:,2), ks);
%%
perf(X,ks);











%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function res = genr(k, pos, neg, prob, K_spa)
[W,N] = size(pos);
if rand < prob
    if ~exist('K_spa', 'var')
        res = randi(2,1,N)-1;
    else
        res = scode(N, K_spa, 1);
    end
else
    if k == 1
        res = pos(randi(W),:);
    else
        res = neg(randi(W),:);
    end
end
end

function code = scode(n, k, w)
code = zeros(w, n);
for i = 1:w
    code(i,randperm(n,k)) = 1;
end
end

function perf(X,ks, ks2, tr_only, l1, l2)
errf = @(k,p) mean(k(:)~=p(:));
if exist('tr_only', 'var') && tr_only
    tr_frac = 1;
else
    tr_frac = 0.7;
end
algs = my_algs({'mvnb2', 'ecoclin'}, {'original', 'shuf'}, true, 0);
PARLOOPS = 4;
train_err = cell(2,numel(algs)); test_err = cell(2,numel(algs));
for i = 1:numel(algs)
    [train_err{1,i}, test_err{1,i}] = evaluate_alg(algs(i),...
        X, ks, 'eval_f', errf,...
        'train_frac', tr_frac, 'par_loops', PARLOOPS);
    if exist('ks2', 'var')
        [train_err{2,i}, test_err{2,i}] = evaluate_alg(algs(i),...
        X, ks2, 'eval_f', errf,...
        'train_frac', tr_frac, 'par_loops', PARLOOPS);
    end
end
if ~exist('ks2', 'var')
    train_err(2,:) = train_err(1,:);
    test_err(2,:) = test_err(1,:);
end
figure;
if exist('l1', 'var') && exist('l2', 'var')
    dayset(1).label = l1;
    dayset(2).label = l2;
else
    dayset(1).label = 'sim';
    dayset(2).label = 'sim2';
end

if exist('tr_only', 'var') && tr_only
    plotmat('error', train_err, train_err, algs, dayset, 1-1/length(unique(ks2)));
    return;
end
plotmat('error', train_err, test_err, algs, dayset, 1-1/length(unique(ks)));
end


function [X, ks] = data_gen(M, pos_words, neg_words)
N = size(pos_words, 2);
ks = randi(2,1,M) - 1;
X = zeros(M,N);
for i = 1:M
    X(i,:) = genr(ks(i), pos_words, neg_words, 0);
end
X_noise = rand(size(X));
X(X_noise > 0.2) = 0; % fires only p of the time
X = sparse(X);
end


%%
