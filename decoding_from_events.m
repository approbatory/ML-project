load('../linear_track/Mouse2022/Mouse-2022-20150326-linear-track/Mouse-2022-20150326_093722-linear-track-TracesAndEvents.mat');
%%
X_full = tracesEvents.rawTraces(91:end,:);
X_spike_full = tracesEvents.spikeDeconv(91:end,:);
X_bin_full = Utils.event_detection(X_full);
y_full = tracesEvents.position(91:end,1);

%%
figure;
n_frames = size(X_full,1);
t_time = (1:n_frames)/20; %in s
c_ix = 100;
plot(t_time, X_full(:,c_ix)); hold on;
plot(t_time, X_bin_full(:,c_ix));


%%
%padding of 2s
%X_bin_full_padded = conv2(X_bin_full, ones(2*20,1), 'full');
%X_bin_full_padded(1:size(X_bin_full,1),:);
X_bin_full_padded = pad_X(X_bin_full, 2);
padding_func = @(X_tr) cellfun(@pad_X, X_tr);
%%
res = decoding_proc(X_bin_full_padded, y_full);
%%
figure;
me_avg_tr = mean(res.me_whole_tr,1);
me_s_avg_tr = mean(res.me_s_whole_tr,1);
me_d_avg_tr = mean(res.me_d_whole_tr,1);
funcs = {@mean, @(x)std(x)./sqrt(size(x,1))};
%funcs = {@mean, @std};
%h_us = shadedErrorBar([], res.me_whole_tr, funcs, 'lineprops', 'b'); hold on;
%h_s = shadedErrorBar([], res.me_s_whole_tr, funcs, 'lineprops', 'r');
%h_d = shadedErrorBar([], res.me_d_whole_tr, funcs, 'lineprops', 'm');
h_us = plot(res.me_whole_tr, 'b'); hold on
h_s = plot(res.me_s_whole_tr, 'r');
h_d = plot(res.me_d_whole_tr, 'm');

l_ = refline(0, median(me_avg_tr)); l_.Color = 'b';
l_ = refline(0, median(me_s_avg_tr)); l_.Color = 'r';
l_ = refline(0, median(me_d_avg_tr)); l_.Color = 'm';
%legend([h_us.mainLine, h_s.mainLine, h_d.mainLine], 'Unshuffled', 'Shuffled', 'Diagonal')
legend([h_us, h_s, h_d], 'Unshuffled', 'Shuffled', 'Diagonal')
set(gca, 'YScale', 'log');
box on;
xlabel 'Trial index'
ylabel 'Mean trial error (cm)'
title(sprintf('Mouse2022: FST events, trial-level decoding,\nshowing median reflines'));

%% Using DecodeTensor:
decode_obj = DecodeTensor(4, 'FST_events');
[me, mse, ps, ks, model] = decode_obj.basic_decode(false, [], []);
[me_s, mse_s, ps_s, ks_s, model_s] = decode_obj.basic_decode(true, [], []);
%%
decode_obj = DecodeTensor(4, 'FST_padded');
[p_me, p_mse, p_ps, p_ks, p_model] = decode_obj.basic_decode(false, [], []);
[p_me_s, p_mse_s, p_ps_s, p_ks_s, p_model_s] = decode_obj.basic_decode(true, [], []);
[p_me_d, p_mse_d, p_ps_d, p_ks_d, p_model_d] = decode_obj.basic_decode(false, [], [], my_algs('ecoclin', 'shuf'));
%%
p_me_tr = mean(abs(reshape(ceil(p_ks/2), 20, []) - reshape(ceil(p_ps/2), 20, []))*5.9);
p_me_s_tr = mean(abs(reshape(ceil(p_ks_s/2), 20, []) - reshape(ceil(p_ps_s/2), 20, []))*5.9);
p_me_d_tr = mean(abs(reshape(ceil(p_ks_d/2), 20, []) - reshape(ceil(p_ps_d/2), 20, []))*5.9);
figure; 
plot(p_me_tr, 'b'); hold on; 
plot(p_me_s_tr, 'r');
plot(p_me_d_tr, 'm');
l_ = refline(0, median(p_me_tr)); l_.Color = 'b';
l_ = refline(0, median(p_me_s_tr)); l_.Color = 'r';
l_ = refline(0, median(p_me_d_tr)); l_.Color = 'm';
legend Unshuffled Shuffled Diagonal
set(gca, 'YScale', 'log');
box on;
xlabel 'Trial index'
ylabel 'Mean trial error (cm)'
title(sprintf('Mouse2022: FST padded (2s), trial-level decoding,\nshowing median reflines'));
%%
function res = decoding_proc(X_full, y_full)
opt = DecodeTensor.default_opt;
[~,~,tr_start,tr_end,~,~,ks] = DecodeTensor.new_sel(y_full, opt);
[X_tr, ks_tr] = split_trials(X_full, ks, tr_start, tr_end);
num_trials = numel(tr_start);
alg = my_algs('ecoclin');
n_reps = 1;
for rep = 1:n_reps
    train_trials = mod(randperm(num_trials),2) == 1;
    test_trials = ~train_trials;
    decode_timing = tic;
    for q = 1:2
        X_tr_train = X_tr(train_trials);
        ks_tr_train = ks_tr(train_trials);
        X_s_tr_train = shuffle_cells(X_tr_train, ks_tr_train);
        model = alg.train(cell2mat(X_tr_train), cell2mat(ks_tr_train));
        model_s = alg.train(cell2mat(X_s_tr_train), cell2mat(ks_tr_train));
        
        
        X_tr_test= X_tr(test_trials);
        ks_tr_test = ks_tr(test_trials);
        X_s_tr_test = shuffle_cells(X_tr_test, ks_tr_test);
        
        used_lengths = cellfun(@(x) size(x,1), ks_tr_test);
        ps_tr = mat2cell(alg.test(model, cell2mat(X_tr_test)), used_lengths);
        ps_s_tr = mat2cell(alg.test(model_s, cell2mat(X_s_tr_test)), used_lengths);
        ps_d_tr = mat2cell(alg.test(model_s, cell2mat(X_tr_test)), used_lengths);
        % evaluation
        me_func = @(k,p) opt.bin_width.*mean(abs(ceil(k/2) - ceil(p/2)));
        mse_func = @(k,p) opt.bin_width.^2.*mean((ceil(k/2) - ceil(p/2)).^2);
        mse_tr = cellfun(mse_func, ks_tr_test, ps_tr);
        mse_s_tr = cellfun(mse_func, ks_tr_test, ps_s_tr);
        mse_d_tr = cellfun(mse_func, ks_tr_test, ps_d_tr);
        
        me_tr = cellfun(me_func, ks_tr_test, ps_tr);
        me_s_tr = cellfun(me_func, ks_tr_test, ps_s_tr);
        me_d_tr = cellfun(me_func, ks_tr_test, ps_d_tr);
        

        mse_whole_tr{rep}(test_trials) = mse_tr;
        mse_s_whole_tr{rep}(test_trials) = mse_s_tr;
        mse_d_whole_tr{rep}(test_trials) = mse_d_tr;
        
        me_whole_tr{rep}(test_trials) = me_tr;
        me_s_whole_tr{rep}(test_trials) = me_s_tr;
        me_d_whole_tr{rep}(test_trials) = me_d_tr;
        
        tmp = train_trials;
        train_trials = test_trials;
        test_trials = tmp;
    end
    n_seconds = toc(decode_timing);
    fprintf('Done rep %d/%d, in time %.2f s\n', rep, n_reps, n_seconds);
end

%
res.mse_whole_tr = cell2mat(mse_whole_tr');
res.mse_s_whole_tr = cell2mat(mse_s_whole_tr');
res.mse_d_whole_tr = cell2mat(mse_d_whole_tr');

res.me_whole_tr = cell2mat(me_whole_tr');
res.me_s_whole_tr = cell2mat(me_s_whole_tr');
res.me_d_whole_tr = cell2mat(me_d_whole_tr');
end



%% util funcs
function [X_tr, ks_tr] = split_trials(X_full, ks, tr_start, tr_end)

num_trials = numel(tr_start);

tr_mask = zeros(numel(ks),1);
tr_lens = zeros(1,numel(tr_start));
X_tr = cell(num_trials, 1);
ks_tr = cell(num_trials, 1);
for tr_i = 1:numel(tr_start)
    s = tr_start(tr_i); e = tr_end(tr_i);
    tr_mask(s:e) = tr_i;
    tr_lens(tr_i) = e - s + 1;
    
    X_tr{tr_i} = X_full(s:e,:);
    ks_tr{tr_i} = ks(s:e);
end
end


function X_s_tr = shuffle_cells(X_tr, ks_tr)
used_lengths = cellfun(@(x) size(x,1), ks_tr);
assert(isequal(used_lengths, cellfun(@(x) size(x,1), X_tr)));

X = cell2mat(X_tr);
ks = cell2mat(ks_tr);
X_s = shuffle(X, ks);

X_s_tr = mat2cell(X_s, used_lengths, size(X_s,2));
end


function X_bin_full_padded = pad_X(X_bin_full, pad_seconds)
X_bin_full_padded = conv2(X_bin_full, ones(pad_seconds*20,1), 'full');
X_bin_full_padded = X_bin_full_padded(1:size(X_bin_full,1),:);
end