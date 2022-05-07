function f2b
%% Figure 2b
% Raw numbers for:
% - LineX
% - LineRealY, LineShuffledY, LinePositionY
rng(0);

data = decode_demo;
make_xlsx(data, 'f2b');
end

%% helper functions
function data = decode_demo
data_source = SessManager.load_special('Mouse2022');
opt = DecodeTensor.default_opt;
load(data_source.source_path);
traces = tracesEvents.events_transients;

[~, ~, trial_start, trial_end, trial_dir, ~, ks] =...
    DecodeTensor.new_sel(tracesEvents.position, opt);
within_trial = zeros(size(ks));
trial_start = trial_start(trial_dir==1);
trial_end = trial_end(trial_dir==1);
n_trials = numel(trial_start);
for i = 1:n_trials
    within_trial(trial_start(i):trial_end(i)) = i;
end

time_coord = (1:numel(ks))./opt.samp_freq;
if true
    first_half = (mod(within_trial,2) == 1) & (within_trial ~= 0);
    second_half = (mod(within_trial,2) == 0) & (within_trial ~= 0);
else
    first_half = (within_trial <= n_trials/2) & (within_trial ~= 0);
    second_half = (within_trial > n_trials/2) & (within_trial ~= 0);
end

ks_first_half = ks(first_half);
ks_second_half = ks(second_half);
X_first_half = traces(first_half, :);
X_second_half = traces(second_half, :);
X_first_half_s = shuffle(X_first_half, ks_first_half);
X_second_half_s = shuffle(X_second_half, ks_second_half);

alg = my_algs('ecoclin');
model1 = alg.train(X_first_half, ks_first_half); disp('m1');
model1_s = alg.train(X_first_half_s, ks_first_half); disp('m1s');
model2 = alg.train(X_second_half, ks_second_half); disp('m2');
model2_s = alg.train(X_second_half_s, ks_second_half); disp('m2s');

ps_first_half = alg.test(model2, X_first_half);
ps_first_half_s = alg.test(model2_s, X_first_half_s);
ps_first_half_d = alg.test(model2_s, X_first_half);
ps_second_half = alg.test(model1, X_second_half);
ps_second_half_s = alg.test(model1_s, X_second_half_s);
ps_second_half_d = alg.test(model1_s, X_second_half);

me_ = @(k,p) mean(abs(ceil(k(:)/2)-ceil(p(:)/2))).*opt.bin_width;
mse_ = @(k,p) mean((ceil(k(:)/2)-ceil(p(:)/2)).^2).*opt.bin_width.^2;
fprintf('Mean error: unsh: %f, sh: %f\n',  me_(ks_first_half, ps_first_half),...
    me_(ks_first_half, ps_first_half_s));
fprintf('Mean error: unsh: %f, sh: %f\n',  me_(ks_second_half, ps_second_half),...
    me_(ks_second_half, ps_second_half_s));


tr_track = within_trial(first_half|second_half);

ks_ = ks(within_trial~=0);
ps_(first_half) = ps_first_half;
ps_(second_half) = ps_second_half;
ps_ = ps_(within_trial~=0)';
ps_s_(first_half) = ps_first_half_s;
ps_s_(second_half) = ps_second_half_s;
ps_s_ = ps_s_(within_trial~=0)';

for i = 1:max(tr_track)
    tr_err(i) = me_(ks_(tr_track==i), ps_(tr_track==i));
    tr_err_s(i) = me_(ks_(tr_track==i), ps_s_(tr_track==i));
end

figure;
t = (1:numel(ks_first_half))./opt.samp_freq;
t_start = 83.5 + 2.05; %which to show
t_end = 83.5 + 4.05;%90;

hold on;
h(1) = plot(t - t_start, (ceil(ps_first_half/2) - 0.5)*opt.bin_width, '-b');
h(2) = plot(t - t_start, (ceil(ps_first_half_s/2) - 0.5)*opt.bin_width, '-r');
h(3) = plot(t - t_start, (ceil(ks_first_half/2) - 0.5)*opt.bin_width, '-k');

data.LineX = t - t_start;
data.LineRealY = (ceil(ps_first_half/2) - 0.5)*opt.bin_width;
data.LineShuffledY = (ceil(ps_first_half_s/2) - 0.5)*opt.bin_width;
data.LinePositionY = (ceil(ks_first_half/2) - 0.5)*opt.bin_width;

trial_boundaries = (find(diff(within_trial(first_half))>0)+0.5)./opt.samp_freq - t_start;
for i = 1:numel(trial_boundaries)
    x_ = trial_boundaries(i);
    line([x_ x_], ylim, 'Color', 'k', 'LineStyle', '--');
end
xlim([t_start t_end] - t_start);

limits = xlim;
data_filter = (data.LineX >= limits(1)) & (data.LineX <= limits(2));
data.LineX = data.LineX(data_filter);
data.LineRealY = data.LineRealY(data_filter);
data.LineShuffledY = data.LineShuffledY(data_filter);
data.LinePositionY = data.LinePositionY(data_filter);

ylim([-Inf Inf]);
xlabel 'Time (s)';
ylabel 'Position (cm)';
legend(h, 'Real', 'Shuffled', 'Target');
legend boxoff
legend location best
figure_format('boxsize', [1 0.7]*1.5);% box on;
end