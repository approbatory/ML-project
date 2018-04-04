L = 16; N = 1000;
points = linspace(0,L,1000);
means = L*rand(1,N);
stds = ones(1,N);
base_lambda = 0.03;
Lambdas = @(point) base_lambda .* exp(-(point - means).^2./(2*stds)); %(position, neuron)


times = 1:100000; speed = 0.01;
trajectory = L*sin(speed*times).^2;
rates = Lambdas(trajectory');

X = double(sparse(poissrnd(rates) > 0));
ks = ceil(trajectory);
errf = @(k,p) mean(abs(k(:)-p(:)));

%%
algs = my_algs({'mvnb2', 'ecoclin'}, {'original', 'shuf'});
PARLOOPS = 4;
train_err = cell(2,numel(algs)); test_err = cell(2,numel(algs));
for i = 1:numel(algs)
    [train_err{1,i}, test_err{1,i}] = evaluate_alg(algs(i),...
        X, ks, 'eval_f', errf,...
        'train_frac', 0.7, 'par_loops', PARLOOPS);
end
train_err(2,:) = train_err(1,:);
test_err(2,:) = test_err(1,:);
figure;
dayset(1).label = 'sim';
dayset(2).label = 'sim2';
%%
plotmat('error', train_err, test_err, algs, dayset, 5);


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
N = 1000;
W = 500;
sparsity = 0.03;
K = floor(N*sparsity);
pos_words = scode(N,K,W);%randi(2,W,N)-1;
neg_words = scode(N,K,W);%randi(2,W,N)-1;
M = 10000;
ks = randi(2,1,M) - 1;
X = zeros(M,N);
for i = 1:M
    X(i,:) = genr(ks(i), pos_words, neg_words, 0.1, K);
end
X = sparse(X);
%X = X .* (rand(size(X)) < 0.1);
errf = @(k,p) mean(k(:)~=p(:));
%%
algs = my_algs({'mvnb2', 'ecoclin'}, {'original', 'shuf'});
PARLOOPS = 64;
train_err = cell(2,numel(algs)); test_err = cell(2,numel(algs));
for i = 1:numel(algs)
    [train_err{1,i}, test_err{1,i}] = evaluate_alg(algs(i),...
        X, ks, 'eval_f', errf,...
        'train_frac', 0.7, 'par_loops', PARLOOPS);
end
train_err(2,:) = train_err(1,:);
test_err(2,:) = test_err(1,:);

figure;
dayset(1).label = 'sim';
dayset(2).label = 'sim2';
plotmat('error', train_err, test_err, algs, dayset, 0.5);












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