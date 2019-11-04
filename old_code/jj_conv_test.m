%try jjconv on linear track:

%loading data
save_file = '../linear_track/pablo_data_ds.mat';

%load the data into a ds array
loaded_res = load(save_file);
my_ds = loaded_res.my_ds;

%%
%arbitrarily choose ix = 7
ix = 7;
ds = my_ds(ix);

alg = my_algs('ecoclin');
num_samples = 64;
midLen = 0.8*120;
num_bins = 20;

using_binary = true;
if using_binary
    my_X = one_time_events(ds.trials.traces).';
else
    my_X = ds.trials.traces.';
end
my_y = ds.trials.centroids;
[sel_fw, sel_bw] = select_directions(my_y);
%just take forwards
my_X = my_X(sel_fw,:);
my_y = my_y(sel_fw,:);

my_binner = @(y) gen_place_bins(y, num_bins, midLen);

my_ticker = tic;
[tr_err, te_err] = evala(alg, my_X, my_y, my_binner,...
    'split', 'nonlocal',...
    'repeats', num_samples, 'verbose', true,...
    'use_par', true, 'errfunc', 'mean_dist');
toc(my_ticker);

%% testing out the convolution

conv_func = @(tau,t) exp(-t./tau).';

%setting tau=10 for now (0.5s)
taus = [0.1 0.2 0.5 1 2 5 10 20 50 100];
for i = 1:numel(taus)
    if using_binary
        X = one_time_events(ds.trials.traces).';
    else
        X = ds.trials.traces.';
    end
    jj_traces = conv2(X, conv_func(taus(i), 0:1000));
    jj_traces = jj_traces(1:size(X,1),:);
    
    my_jj_traces = jj_traces(sel_fw,:);
    fprintf('Using tau = %f :\n', taus(i));
    my_ticker = tic;
    [tr_err_jj{i}, te_err_jj{i}] = evala(alg, my_jj_traces, my_y, my_binner,...
        'split', 'nonlocal',...
        'repeats', num_samples, 'verbose', true,...
        'use_par', true, 'errfunc', 'mean_dist');
    toc(my_ticker);
end

%% plotting
figure;
errorbar(taus/20, cellfun(@mean, te_err_jj), cellfun(@std, te_err_jj)./sqrt(num_samples));
set(gca, 'XScale', 'log');
title('Effect of exponential decay tau on decoding performance');
xlabel('tau (seconds)');
ylabel('Place decoding mean test error (cm)');