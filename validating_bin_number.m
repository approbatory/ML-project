function validating_bin_number(random_str)
my_rand_seed = str2double(random_str);
rng(my_rand_seed);
fprintf('random seed %d\n\n', my_rand_seed);
total_ticker = tic;
load('../linear_track/Mouse2022/Mouse-2022-20150326-linear-track/Mouse-2022-20150326_093722-linear-track-TracesAndEvents.mat');
%%
X_full = tracesEvents.rawTraces(91:end,:);
%X_spike_full = tracesEvents.spikeDeconv(91:end,:);
%X_bin_full = Utils.event_detection(X_full);
%X_bin_full_padded = conv2(X_bin_full, ones(0.8*20,1), 'full');
%X_bin_full_padded = X_bin_full_padded(1:size(X_full,1),:);
%X_full = X_bin_full_padded;
y_full = tracesEvents.position(91:end,1);
%%
nb_list = [5 10 15 18 20 22 25 30 35 40];
for nb_i = 1:numel(nb_list)
opt = DecodeTensor.default_opt;
opt.n_bins = nb_list(nb_i); %CHANGE HERE IN LOOP
opt.bin_width = 118/opt.n_bins;
opt.discard_incomplete_trials = false;
[cpp,~,tr_start,tr_end,~,~,ks] = DecodeTensor.new_sel(y_full, opt);
num_trials = numel(tr_start);

tr_mask = zeros(numel(ks),1);
tr_lens = zeros(1,numel(tr_start));
X_tr = cell(num_trials, 1);
ks_tr = cell(num_trials, 1);
ys_tr = cell(num_trials, 1);
for tr_i = 1:numel(tr_start)
    s = tr_start(tr_i); e = tr_end(tr_i);
    tr_mask(s:e) = tr_i;
    tr_lens(tr_i) = e - s + 1;
    
    X_tr{tr_i} = X_full(s:e,:);
    ks_tr{tr_i} = ks(s:e);
    ys_tr{tr_i} = y_full(s:e);
end
%%
n_reps = 1;
%mse_whole_tr = zeros(n_reps, num_trials);
%mse_s_whole_tr = zeros(n_reps, num_trials);
%mse_d_whole_tr = zeros(n_reps, num_trials);
alg = my_algs('ecoclin');
bin_center_coords = zeros(2*opt.n_bins, 1);
for b_i = 1:2*opt.n_bins
    bin_center_coords(b_i) = mean(y_full(ks == b_i));
end
clear mse_whole_tr mse_s_whole_tr mse_d_whole_tr me_whole_tr me_s_whole_tr me_d_whole_tr
for rep = 1:n_reps
    train_trials = mod(randperm(num_trials),2) == 1;
    test_trials = ~train_trials;
    decode_timing = tic;
    for q = 1:2
        X_tr_train = X_tr(train_trials);
        ks_tr_train = ks_tr(train_trials);
        %X_s_tr_train = shuffle_cells(X_tr_train, ks_tr_train);
        model = alg.train(cell2mat(X_tr_train), cell2mat(ks_tr_train));
        %model_s = alg.train(cell2mat(X_s_tr_train), cell2mat(ks_tr_train));
        
        
        X_tr_test= X_tr(test_trials);
        ks_tr_test = ks_tr(test_trials);
        ys_tr_test = ys_tr(test_trials);
        %X_s_tr_test = shuffle_cells(X_tr_test, ks_tr_test);
        
        used_lengths = cellfun(@(x) size(x,1), ks_tr_test);
        ps_tr = mat2cell(alg.test(model, cell2mat(X_tr_test)), used_lengths);
        %ps_s_tr = mat2cell(alg.test(model_s, cell2mat(X_s_tr_test)), used_lengths);
        %ps_d_tr = mat2cell(alg.test(model_s, cell2mat(X_tr_test)), used_lengths);
        % evaluation
        me_func = @(y,p) cpp.*mean(abs(y - bin_center_coords(p)));
        mse_func = @(y,p) cpp.^2.*mean((y - bin_center_coords(p)).^2);
        mse_tr = cellfun(mse_func, ys_tr_test, ps_tr);
        %mse_s_tr = cellfun(mse_func, ks_tr_test, ps_s_tr);
        %mse_d_tr = cellfun(mse_func, ks_tr_test, ps_d_tr);
        
        me_tr = cellfun(me_func, ys_tr_test, ps_tr);
        %me_s_tr = cellfun(me_func, ks_tr_test, ps_s_tr);
        %me_d_tr = cellfun(me_func, ks_tr_test, ps_d_tr);
        
        %% TODO RUN THIS THING
        
        %mse_whole_tr(rep, test_trials) = mse_tr;
        %mse_s_whole_tr(rep, test_trials) = mse_s_tr;
        %mse_d_whole_tr(rep, test_trials) = mse_d_tr;
        mse_whole_tr{rep}(test_trials) = mse_tr;
        %mse_s_whole_tr{rep}(test_trials) = mse_s_tr;
        %mse_d_whole_tr{rep}(test_trials) = mse_d_tr;
        
        me_whole_tr{rep}(test_trials) = me_tr;
        %me_s_whole_tr{rep}(test_trials) = me_s_tr;
        %me_d_whole_tr{rep}(test_trials) = me_d_tr;
        
        tmp = train_trials;
        train_trials = test_trials;
        test_trials = tmp;
    end
    n_seconds = toc(decode_timing);
    fprintf('Done rep %d/%d, in time %.2f s\n', rep, n_reps, n_seconds);
end

mse_whole_tr = cell2mat(mse_whole_tr');
%mse_s_whole_tr = cell2mat(mse_s_whole_tr');
%mse_d_whole_tr = cell2mat(mse_d_whole_tr');

me_whole_tr = cell2mat(me_whole_tr');
%me_s_whole_tr = cell2mat(me_s_whole_tr');
%me_d_whole_tr = cell2mat(me_d_whole_tr');

me_mean_tr(nb_i) = mean(me_whole_tr);
mse_mean_tr(nb_i) = mean(mse_whole_tr);
end

save(sprintf('validating_bin_number_record_%d.mat', randi(1e9)), 'nb_list', 'me_mean_tr', 'mse_mean_tr');
toc(total_ticker);
exit;