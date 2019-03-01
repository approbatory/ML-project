%non averaged decoding
o = Analyzer('../linear_track/Mouse2022/Mouse-2022-20150326-linear-track/Mouse-2022-20150326_093722-linear-track-TracesAndEvents.mat');
opt = DecodeTensor.default_opt;
[~,~,tr_start,tr_end,~,~,ks] = DecodeTensor.new_sel(o.data.y.raw.full, opt);
num_trials = numel(tr_start);

tr_mask = zeros(numel(ks),1);
tr_lens = zeros(1,numel(tr_start));
X_tr = cell(num_trials, 1);
ks_tr = cell(num_trials, 1);
for tr_i = 1:numel(tr_start)
    s = tr_start(tr_i); e = tr_end(tr_i); 
    tr_mask(s:e) = tr_i;
    tr_lens(tr_i) = e - s + 1;
    
    X_tr{tr_i} = o.data.X.full(s:e,:);
    ks_tr{tr_i} = ks(s:e);
end

train_trials = mod(randperm(num_trials),2) == 1;
test_trials = ~train_trials;

alg = my_algs('ecoclin');

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
%% evaluation
me_tr = cellfun(@(k,p) mean(abs(ceil(k/2) - ceil(p/2))), ks_tr_test, ps_tr);
me_s_tr = cellfun(@(k,p) mean(abs(ceil(k/2) - ceil(p/2))), ks_tr_test, ps_s_tr);

%%
figure;
plot(me_tr*5.9, 'b'); hold on;
plot(me_s_tr*5.9, 'r');
l_ = refline(0, median(me_tr*5.9)); l_.Color = 'b';
l_ = refline(0, median(me_s_tr*5.9)); l_.Color = 'r';
%%
figure;
t_time = ((1:numel(cell2mat(ks_tr_test))) - 1000)/20; %20Hz
plot(t_time, (ceil(cell2mat(ks_tr_test)/2)-0.5)*5.9, 'k'); hold on;
plot(t_time, (ceil(cell2mat(ps_tr)/2)-0.5)*5.9, 'b');
plot(t_time, (ceil(cell2mat(ps_s_tr)/2)-0.5)*5.9, 'r');
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