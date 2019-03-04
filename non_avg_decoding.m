load('../linear_track/Mouse2022/Mouse-2022-20150326-linear-track/Mouse-2022-20150326_093722-linear-track-TracesAndEvents.mat');
%%
X_full = tracesEvents.rawTraces(91:end,:);
X_spike_full = tracesEvents.spikeDeconv(91:end,:);
X_bin_full = Utils.event_detection(X_full);
y_full = tracesEvents.position(91:end,1);
%%
opt = DecodeTensor.default_opt;
[~,~,tr_start,tr_end,~,~,ks] = DecodeTensor.new_sel(y_full, opt);
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

n_reps = 20;
%mse_whole_tr = zeros(n_reps, num_trials);
%mse_s_whole_tr = zeros(n_reps, num_trials);
%mse_d_whole_tr = zeros(n_reps, num_trials);
alg = my_algs('ecoclin');
parfor rep = 1:n_reps
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
        
        
        
        %mse_whole_tr(rep, test_trials) = mse_tr;
        %mse_s_whole_tr(rep, test_trials) = mse_s_tr;
        %mse_d_whole_tr(rep, test_trials) = mse_d_tr;
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
mse_whole_tr = cell2mat(mse_whole_tr');
mse_s_whole_tr = cell2mat(mse_s_whole_tr');
mse_d_whole_tr = cell2mat(mse_d_whole_tr');

me_whole_tr = cell2mat(me_whole_tr');
me_s_whole_tr = cell2mat(me_s_whole_tr');
me_d_whole_tr = cell2mat(me_d_whole_tr');
return;
%%
figure;
me_avg_tr = mean(me_whole_tr);
me_s_avg_tr = mean(me_s_whole_tr);
me_d_avg_tr = mean(me_d_whole_tr);
funcs = {@mean, @(x)std(x)./sqrt(size(x,1))};
%funcs = {@mean, @std};
h_us = shadedErrorBar([], me_whole_tr, funcs, 'lineprops', 'b'); hold on;
h_s = shadedErrorBar([], me_s_whole_tr, funcs, 'lineprops', 'r');
h_d = shadedErrorBar([], me_d_whole_tr, funcs, 'lineprops', 'm');
l_ = refline(0, median(me_avg_tr)); l_.Color = 'b';
l_ = refline(0, median(me_s_avg_tr)); l_.Color = 'r';
l_ = refline(0, median(me_d_avg_tr)); l_.Color = 'm';
legend([h_us.mainLine, h_s.mainLine, h_d.mainLine], 'Unshuffled', 'Shuffled', 'Diagonal')
set(gca, 'YScale', 'log');
box on;
xlabel 'Trial index'
ylabel 'Mean trial error (cm)'
title(sprintf('Mouse2022: rawTraces, trial-level decoding,\nshowing median reflines'));
%%
figure;
t_time = ((1:numel(cell2mat(ks_tr_test))) - 1000)/20; %20Hz
plot(t_time, (ceil(cell2mat(ks_tr_test)/2)-0.5)*5.9, 'k'); hold on;
plot(t_time, (ceil(cell2mat(ps_tr)/2)-0.5)*5.9, 'b');
plot(t_time, (ceil(cell2mat(ps_s_tr)/2)-0.5)*5.9, 'r');
plot(t_time, (ceil(cell2mat(ps_d_tr)/2)-0.5)*5.9, 'm');
trial_marks = (cumsum(cellfun(@(x)size(x,1), ks_tr_test)) - 1000)/20; %20Hz
for i = 1:numel(trial_marks)
    m = trial_marks(i);
    line([m m], ylim, 'Color', 'black', 'LineStyle', '--');
end
xlim([0 9]);
ylim([0 118]);
xlabel 'Time (s)'
ylabel 'Position (cm)'

figure_format('boxsize', [0.8 0.7]*1.05); box on;
%%
function X_s_tr = shuffle_cells(X_tr, ks_tr)
used_lengths = cellfun(@(x) size(x,1), ks_tr);
assert(isequal(used_lengths, cellfun(@(x) size(x,1), X_tr)));

X = cell2mat(X_tr);
ks = cell2mat(ks_tr);
X_s = shuffle(X, ks);

X_s_tr = mat2cell(X_s, used_lengths, size(X_s,2));
end