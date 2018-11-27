function test_selector(dispatch_index)
source_paths = {'../linear_track/Mouse2010/Mouse-2010-20141125-linear-track/Mouse-2010-20141125_000156-linear-track-TracesAndEvents.mat',...
    '../linear_track/Mouse2012/Mouse-2012-20150114-linear-track/Mouse-2012-20150114_140933-linear-track-TracesAndEvents.mat',...
    '../linear_track/Mouse2019/Mouse-2019-20150311-linear-track/Mouse-2019-20150311_101049-linear-track-TracesAndEvents.mat',...
    '../linear_track/Mouse2022/Mouse-2022-20150326-linear-track/Mouse-2022-20150326_093722-linear-track-TracesAndEvents.mat',...
    '../linear_track/Mouse2023/Mouse-2023-20150326-linear-track/Mouse-2023-20150326_101329-linear-track-TracesAndEvents.mat',...
    '../linear_track/Mouse2024/Mouse-2024-20150311-linear-track/Mouse-2024-20150311_073912-linear-track-TracesAndEvents.mat',...
    '../linear_track/Mouse2026/Mouse-2026-20150303-linear-track/Mouse-2026-20150303_095846-linear-track-TracesAndEvents.mat',...
    '../linear_track/Mouse2028/Mouse-2028-20150327-linear-track/Mouse-2028-20150327_105544-linear-track-TracesAndEvents.mat'};

mouse_names = {'Mouse2010', 'Mouse2012', 'Mouse2019', 'Mouse2022',...
    'Mouse2023', 'Mouse2024', 'Mouse2026', 'Mouse2028'};

opt.total_length = 118; %cm
opt.cutoff_p = 5; %percentile
opt.samp_freq = 20; %Hz
opt.v_thresh = 4; %cm/s
opt.n_bins = 20;
opt.d_neurons = 30;

opt.bin_width = opt.total_length/opt.n_bins;

decode_series(source_paths{dispatch_index}, mouse_names{dispatch_index},...
    opt, 'tensor.db');
end

function decode_series(source_path, mouse_id, opt, database_file)
if ~exist(database_file, 'file')
    conn = sqlite(database_file, 'create');
    conn.exec('create table decoding (Mouse text, Setting text, NumNeurons int, DataSize int, MeanErrors real, MSE real)');
else
    conn = sqlite(database_file);
end

table_name = 'decoding';
field_names = {'Mouse', 'Setting', 'NumNeurons', 'DataSize', 'MeanErrors', 'MSE'};
[data_tensor, tr_dir] = tensor_loader(source_path, mouse_id, opt);

num_neurons = size(data_tensor, 1);
num_trials = min(sum(tr_dir == 1), sum(tr_dir == -1));

alg = my_algs('ecoclin');
alg_diag = my_algs('ecoclin', 'shuf');

neuron_series = opt.d_neurons:opt.d_neurons:num_neurons;
if neuron_series(1) ~= 1
    neuron_series = [1 neuron_series];
end
if neuron_series(end) ~= num_neurons
    neuron_series = [neuron_series num_neurons];
end


for n_neu = neuron_series
    [mean_err, MSE] = decode_tensor(data_tensor, tr_dir, opt.bin_width, alg, false,...
        n_neu, num_trials);
    conn.insert(table_name, field_names,...
        {mouse_id, 'unshuffled', n_neu, num_trials, mean_err, MSE});
    fprintf('n_neu=%d\tmean_err = %.2f\n', n_neu, mean_err);
    
    [mean_err_s, MSE_s] = decode_tensor(data_tensor, tr_dir, opt.bin_width, alg, true,...
        n_neu, num_trials);
    conn.insert(table_name, field_names,...
        {mouse_id, 'shuffled', n_neu, num_trials, mean_err_s, MSE_s});
    fprintf('n_neu=%d\tmean_err_s = %.2f\n', n_neu, mean_err_s);
    
    [mean_err_d, MSE_d] = decode_tensor(data_tensor, tr_dir, opt.bin_width, alg_diag, false,...
        n_neu, num_trials);
    conn.insert(table_name, field_names,...
        {mouse_id, 'diagonal', n_neu, num_trials, mean_err_d, MSE_d});
    fprintf('n_neu=%d\tmean_err_d = %.2f\n\n', n_neu, mean_err_d);
end
conn.close;
end

function [T, d] = tensor_loader(source_path, mouse_id, opt)
load(source_path);
track_coord = tracesEvents.position(:,1);
X = tracesEvents.rawTraces;
if strcmp(mouse_id, 'Mouse2022')
    track_coord = track_coord(91:end);
    X = X(91:end,:);
end

[~, ~, tr_s, tr_e, tr_dir, tr_bins, ~] = new_sel(track_coord, opt);
data_tensor = construct_tensor(X, tr_bins, opt.n_bins, tr_s, tr_e);
T = data_tensor; d = tr_dir;
end

function [mean_err, MSE, ps, ks] = decode_tensor(data_tensor, tr_dir,...
    binsize, alg, shuf, num_neurons, num_trials)
[data_tensor, tr_dir] = cut_tensor(data_tensor, tr_dir, num_neurons, num_trials);
[~, n_bins, ~] = size(data_tensor);
[T1, d1, T2, d2, division] = holdout_half(data_tensor, tr_dir);
if shuf
    T1 = shuffle_tensor(T1, d1);
    T2 = shuffle_tensor(T2, d2);
end
[sup_X1, sup_ks1] = tensor2dataset(T1, d1);
[sup_X2, sup_ks2] = tensor2dataset(T2, d2);

mean_err_func = @(ks, ps) mean(abs(ceil(ks/2) - ceil(ps/2))) * binsize;
MSE_func = @(ks, ps) mean((ceil(ks/2) - ceil(ps/2)).^2) * binsize;

model = alg.train(sup_X1, sup_ks1);
sup_ps2 = alg.test(model, sup_X2);
mean_err2 = mean_err_func(sup_ks2, sup_ps2);
MSE2 = MSE_func(sup_ks2, sup_ps2);

model = alg.train(sup_X2, sup_ks2);
sup_ps1 = alg.test(model, sup_X1);
mean_err1 = mean_err_func(sup_ks1, sup_ps1);
MSE1 = MSE_func(sup_ks1, sup_ps1);

mean_err = mean([mean_err1 mean_err2]);
MSE = mean([MSE1 MSE2]);

ps( repmat(division,1,n_bins)) = sup_ps1;
ps(~repmat(division,1,n_bins)) = sup_ps2;

ks( repmat(division,1,n_bins)) = sup_ks1;
ks(~repmat(division,1,n_bins)) = sup_ks2;

[~, tot_ks] = tensor2dataset(data_tensor, tr_dir);
assert(isequal(ks(:), tot_ks(:)));%%TODO remove
end

function [T_cut, d_cut] = cut_tensor(data_tensor, tr_dir, num_neurons, num_trials)
[tot_neurons, ~, ~] = size(data_tensor);

if isempty(num_neurons)
    num_neurons = tot_neurons;
end
assert(num_neurons <= tot_neurons, 'The size of the subset of neurons must be <= the total recorded neurons');
num_right_trials = sum(tr_dir == 1);
num_left_trials = sum(tr_dir == -1);
if isempty(num_trials)
    num_trials = min(num_right_trials, num_left_trials);
end
assert(num_trials <= min(num_right_trials, num_left_trials),...
    'The size of the subset of trials must be <= the total recorded trials');

neuron_subset = randperm(tot_neurons) <= num_neurons;

right_trials_subset = randperm(num_right_trials) <= num_trials;
left_trials_subset = randperm(num_left_trials) <= num_trials;
total_trials_subset(tr_dir == 1) = right_trials_subset;
total_trials_subset(tr_dir ==-1) = left_trials_subset;

T_cut = data_tensor(neuron_subset,:,total_trials_subset);
d_cut = tr_dir(total_trials_subset);
end

function [T1, d1, T2, d2, division] = holdout_half(data_tensor, tr_dir)
num_right_trials = sum(tr_dir == 1);
num_left_trials = sum(tr_dir == -1);
div_r = randperm(num_right_trials) <= num_right_trials/2;
div_l = randperm(num_left_trials) <= num_left_trials/2;

division(tr_dir == 1) = div_r;
division(tr_dir ==-1) = div_l;

T1 = data_tensor(:,:, division);
T2 = data_tensor(:,:,~division);
d1 = tr_dir( division);
d2 = tr_dir(~division);
end

function [sup_X, sup_ks] = tensor2dataset(data_tensor, tr_dir)
[n_neurons, n_bins, n_trials] = size(data_tensor);
sup_X = zeros(n_bins*n_trials, n_neurons);
sup_ks = zeros(n_bins*n_trials,1);
for t = 1:n_trials
    for b = 1:n_bins
        sup_X((b-1)*n_trials + t,:) = data_tensor(:, b, t);
        sup_ks((b-1)*n_trials + t) = 2*b - (tr_dir(t) == 1);
    end
end
end

function shuf_tensor = shuffle_tensor(data_tensor, tr_dir)
right_tensor = data_tensor(:,:,tr_dir == 1);
num_right_trials = sum(tr_dir == 1);
left_tensor = data_tensor(:,:,tr_dir == -1);
num_left_trials = sum(tr_dir == -1);

[n_neurons, n_bins, ~] = size(data_tensor);
for n = 1:n_neurons
    for b = 1:n_bins
        right_tensor(n,b,:) = right_tensor(n,b, randperm(num_right_trials));
        left_tensor(n,b,:) = left_tensor(n,b, randperm(num_left_trials));
    end
end
shuf_tensor(:,:, tr_dir == 1) = right_tensor;
shuf_tensor(:,:, tr_dir == -1) = left_tensor;
end

function [data_tensor, counts_mat] = construct_tensor(X, tr_bins, n_bins, tr_s, tr_e)
[~, n_neurons] = size(X);
n_trials = length(tr_s);
data_tensor = zeros(n_neurons, n_bins, n_trials);
counts_mat = zeros(n_bins, n_trials);
for tr_i = 1:n_trials
    trial_slice = tr_s(tr_i):tr_e(tr_i);
    trial_data = X(trial_slice, :);
    trial_bins = tr_bins(trial_slice);
    for b = 1:n_bins
        data_tensor(:, b, tr_i) =...
            mean(trial_data(trial_bins == b,:),1);
        counts_mat(b, tr_i) = size(trial_data(trial_bins == b,:),1);
    end
end
end

function [cpp, vel, trial_start, trial_end,...
    trial_direction, track_bins, track_dir_bins] = new_sel(XY, opt)
total_length = opt.total_length; %cm
cutoff_p = opt.cutoff_p; %percentile
samp_freq = opt.samp_freq; %Hz
v_thresh = opt.v_thresh; %cm/s
n_bins = opt.n_bins;
leeway_frac = 1/n_bins;

track_coord = XY(:,1);
track_range = (prctile(track_coord, 100-cutoff_p) -...
    prctile(track_coord, cutoff_p));
cpp = total_length / track_range;
vel = [0; diff(track_coord)] .* cpp .* samp_freq;

fast_frames = abs(vel) > v_thresh;
trial_start = find(diff(fast_frames) == 1);
trial_end = find(diff(fast_frames) == -1);

sc = track_coord(trial_start); ec = track_coord(trial_end);
b = prctile(track_coord, cutoff_p) + track_range*leeway_frac;
t = prctile(track_coord, 100-cutoff_p) - track_range*leeway_frac;
regular_trial = ((sc < b) & (ec > t)) | ((ec < b) & (sc > t));

trial_start = trial_start(regular_trial);
trial_end = trial_end(regular_trial);

trial_direction((track_coord(trial_start) < b) & (track_coord(trial_end) > t)) = 1;
trial_direction((track_coord(trial_end) < b) & (track_coord(trial_start) > t)) = -1;


binner = @(y, n_bins)...
    ceil(n_bins.*(y - prctile(track_coord, cutoff_p))./track_range);

add_in_direction = @(bins, vel) -(sign(vel)==1) + 2.*bins;

track_bins = binner(track_coord, n_bins);

track_bins(track_bins < 1) = 1;
track_bins(track_bins > n_bins) = n_bins;

track_dir_bins = add_in_direction(track_bins, vel);
end