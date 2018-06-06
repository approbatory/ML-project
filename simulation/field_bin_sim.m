function res = field_bin_sim(n_bins, n_fields)
L = 1;
N = 1000;
means = sort(L*rand(N,1), 'descend');
stds = ones(N,1)*L/n_fields;
base_lambda = 0.05;
Lambdas = @(point) base_lambda .* ...
    (exp(-(point - means).^2./(2*stds.^2)));

times = 1:10000;
speed = 0.01;
trajectory = L*sin(speed*times).^2;
trajectory(trajectory == 0) = eps;
ks = ceil(trajectory ./ (L./n_bins));

rates = Lambdas(trajectory);
X = double(sparse(rand(size(rates)) < rates))';

errf = @(k,p) mean(k(:)~=p(:));

algs = my_algs({'ecoclin'}, {'original', 'shuf'}, true, 0);
par_loops = 16;
train_frac = 0.7;
for i = 1:numel(algs)
    [res.train{i}, res.test{i}] = evaluate_alg(algs(i),...
        X, ks, 'eval_f', errf,...
        'train_frac', train_frac, 'par_loops', par_loops, 'verbose', false);
end