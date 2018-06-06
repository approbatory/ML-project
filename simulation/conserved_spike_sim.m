function res = conserved_spike_sim(n_bins, n_fields)
%n_fields = 16;
%n_bins = 16;

L = 100; %cm
N = 1000;
means = sort(L*rand(N,1));
stds = ones(N,1)*L/n_fields;

spikes_per_lap = 10;
total_time = 1000; %seconds
delta_t = 0.05; %seconds
total_frames = total_time/delta_t;
speed = 20; %cm/s
delta_x = speed*delta_t; %cm
frames_per_lap = L/delta_x;
num_laps = total_frames / frames_per_lap;
%%
spikes = zeros(N, num_laps, spikes_per_lap);
for i = 1:N
    pd = makedist('Normal', 'mu', means(i), 'sigma', stds(i));
    pd = truncate(pd, 0, L);
    spikes(i,:,:) = random(pd, num_laps, spikes_per_lap);
end

Xs = cell(num_laps,1);
kss = cell(1,num_laps);
for i = 1:num_laps
    my_spikes = spikes(:,i,:);
    [neuron, ~, place] = find(my_spikes);
    [s_place, ord] = sort(place);
    s_neuron = neuron(ord);
    
    
    time_bins = ceil(s_place./delta_x);
    Xs{i} = double(sparse(time_bins, s_neuron, 1) ~= 0);
    kss{i} = ceil((1:max(time_bins)).*delta_x ./ (L./n_bins));
end

ks = cell2mat(kss);
X = cell2mat(Xs);
%%
errf = @(k,p) mean(k(:)~=p(:));

%algs = my_algs({'ecoclin'}, {'original', 'shuf'},...
%    true, 0);
algs = my_algs({'mvnb2', 'ecoclin'}, {'original', 'shuf'},...
    true, 0);
par_loops = 16;
train_frac = 0.7;
for i = 1:numel(algs)
    [res.train{i}, res.test{i}] = evaluate_alg(algs(i),...
        X, ks, 'eval_f', errf,...
        'train_frac', train_frac, 'par_loops', par_loops, 'verbose', true);
end
return;
%%
algs = my_algs({'mvnb2', 'ecoclin'}, {'original', 'shuf'},...
    true, 0);
LL = 0.01;
algs_reg = my_algs({'mvnb2', 'ecoclin'}, {'original', 'shuf'},...
    true, LL);
par_loops = 16;
train_frac = 0.7;
for i = 1:numel(algs)
    [train_err{1,i}, test_err{1,i}] = evaluate_alg(algs(i),...
        X, ks, 'eval_f', errf,...
        'train_frac', train_frac, 'par_loops', par_loops, 'verbose', true);
    [train_err{2,i}, test_err{2,i}] = evaluate_alg(algs_reg(i),...
        X, ks, 'eval_f', errf,...
        'train_frac', train_frac, 'par_loops', par_loops, 'verbose', true);
    [train_err{3,i}, test_err{3,i}] = evaluate_alg(algs(i),...
        @()shuffle(X,ks), ks, 'eval_f', errf,...
        'train_frac', train_frac, 'par_loops', par_loops, 'X_is_func', true, 'verbose', true);    
end
%%
berr({'no reg.', sprintf('%f L1 reg.',LL), 'preshuf, noreg'}, test_err, train_err, {algs.name});