top_dir = '../fabio_data';

behavior = {'dhpc08_behavior.txt', 'dhpc09_behavior.txt', 'dhpc10_behavior.txt'};
neural_data = {'dhpc08_cnmfe_spikes.txt', 'dhpc09_cnmfe.mat', 'dhpc10_cnmfe.mat'};

f_ = @(x) fullfile(top_dir, x);

for i = 1:numel(neural_data)
    t_ = tic;
    fs_sess{i} = OpenField(f_(behavior{i}),...
        f_(neural_data{i}));
    parfor j = 1:20
        megares_conv(j) = real_soc_cos(fs_sess{i}, true);
        megares_noconv(j) = real_soc_cos(fs_sess{i}, false);
    end
    res_fs(i).conv = megares_conv;
    res_fs(i).noconv = megares_noconv;
    fprintf('Done fs sess %d in %g s\n', i, toc(t_));
end

sess_idx = SessManager.special_sessions_usable_index({'Mouse2022', 'Mouse2024', 'Mouse2028'});

for i = 1:numel(sess_idx)
    t_ = tic;
    parfor j = 1:20
        res_linear_conv(j) = our_results(sess_idx(i), true);
        res_linear_noconv(j) = our_results(sess_idx(i), false);
    end
    res_us(i).conv = res_linear_conv;
    res_us(i).noconv = res_linear_noconv;
    fprintf('Done us sess %d in %g s\n', i, toc(t_));
end

save stef_fig_compute.mat fs_sess res_fs sess_idx res_us


function megares = real_soc_cos(fs_sess, with_conv)
S = full(fs_sess.spike_traces.');
trim_mobile = @(s) s(fs_sess.mobile, :);

if ~with_conv
    Xm = trim_mobile(S);
else
    Xm = trim_mobile(fs_sess.convolve_data(S));
end

[~, ks_raw] = fs_sess.discrete_pos;

ksm = -ks_raw.*(-1).^fs_sess.mobile;
if ~with_conv
    Xm_cos  = @() trim_mobile(fshuffle(S, ksm, false));
    Xm_coas = @() trim_mobile(fshuffle(S, ksm, true ));
    
    Xm_soc  = @() trim_mobile(fshuffle(S, ksm, false));
    Xm_asoc = @() trim_mobile(fshuffle(S, ksm, true ));
else
    Xm_cos  = @() trim_mobile(fs_sess.convolve_data(fshuffle(S, ksm, false)));
    Xm_coas = @() trim_mobile(fs_sess.convolve_data(fshuffle(S, ksm, true )));
    
    Xm_soc  = @() trim_mobile(fshuffle(fs_sess.convolve_data(S), ksm, false));
    Xm_asoc = @() trim_mobile(fshuffle(fs_sess.convolve_data(S), ksm, true ));
end

fprintf('Getting data...');
[~, y_cont, ks] = fs_sess.get_continuous_dataset;

fprintf('Decoding real...');
t_ = tic;
res = OpenField.decode_all_general(Xm, y_cont, ks, fs_sess.bin_nums, fs_sess.box_dims);
megares.real = res.fs_metric_real;
fprintf(' Done in %g s\n', toc(t_));

fprintf('Decoding cos...');
t_ = tic;
res = OpenField.decode_all_general(Xm_cos(), y_cont, ks, fs_sess.bin_nums, fs_sess.box_dims);
megares.cos = res.fs_metric_real;
fprintf(' Done in %g s\n', toc(t_));

fprintf('Decoding coas...');
t_ = tic;
res = OpenField.decode_all_general(Xm_coas(), y_cont, ks, fs_sess.bin_nums, fs_sess.box_dims);
megares.coas = res.fs_metric_real;
fprintf(' Done in %g s\n', toc(t_));

fprintf('Decoding soc...');
t_ = tic;
res = OpenField.decode_all_general(Xm_soc(), y_cont, ks, fs_sess.bin_nums, fs_sess.box_dims);
megares.soc = res.fs_metric_real;
fprintf(' Done in %g s\n', toc(t_));

fprintf('Decoding asoc...');
t_ = tic;
res = OpenField.decode_all_general(Xm_asoc(), y_cont, ks, fs_sess.bin_nums, fs_sess.box_dims);
megares.asoc = res.fs_metric_real;
fprintf(' Done in %g s\n', toc(t_));

fprintf('Decoding labelshuf...');
t_ = tic;
res = OpenField.decode_all_general(Xm, y_cont, ks(randperm(numel(ks))), fs_sess.bin_nums, fs_sess.box_dims);
megares.labelshuf = res.fs_metric_real;
fprintf(' Done in %g s\n', toc(t_));
end


function res = our_results(usable_idx, with_conv)
bin_width = 5.9;
tau = 0.5; %2.2
fps = 20;
sm = SessManager;
d = sm.cons_usable(usable_idx);
load(d.source_path);
S = tracesEvents.HD_spikes;
if with_conv
    C = tau_conv(S, tau, fps);
else
    C = S;
end

opt = DecodeTensor.default_opt;
[~, ~, tr_s, tr_e, tr_dir, tr_bins, tr_dir_bins] = DecodeTensor.new_sel(tracesEvents.position(:,1), opt);
data_tensor = DecodeTensor.construct_tensor(C, tr_bins, opt.n_bins, tr_s, tr_e);
d.data_tensor = data_tensor;
d.tr_dir = tr_dir;
d.opt = opt;

err_res = DecodeTensor.decode_all(d.data_tensor, d.tr_dir, bin_width, my_algs('ecoclin'), [], []);
res.real = err_res.mean_err.unshuffled;
res.inshuf = err_res.mean_err.shuffled;

within_trial = false(size(tr_dir_bins));
for i = 1:numel(tr_s)
    within_trial(tr_s(i):tr_e(i)) = true;
end

tr_dir_bins_intrial = tr_dir_bins .* (-1).^within_trial;
S_shuf = fshuffle(S, tr_dir_bins_intrial, false);
if with_conv
    C_shuf = tau_conv(S_shuf, tau, fps);
else
    C_shuf = S_shuf;
end
data_tensor = DecodeTensor.construct_tensor(C_shuf, tr_bins, opt.n_bins, tr_s, tr_e);
d.data_tensor = data_tensor;
d.tr_dir = tr_dir;
d.opt = opt;

err_res = DecodeTensor.decode_all(d.data_tensor, d.tr_dir, bin_width, my_algs('ecoclin'), [], []);
res.cos = err_res.mean_err.unshuffled;

tr_dir_bins_intrial = tr_dir_bins .* (-1).^within_trial;
S_shuf = fshuffle(S, tr_dir_bins_intrial, true);
if with_conv
    C_shuf = tau_conv(S_shuf, tau, fps);
else
    C_shuf = S_shuf;
end
data_tensor = DecodeTensor.construct_tensor(C_shuf, tr_bins, opt.n_bins, tr_s, tr_e);
d.data_tensor = data_tensor;
d.tr_dir = tr_dir;
d.opt = opt;

err_res = DecodeTensor.decode_all(d.data_tensor, d.tr_dir, bin_width, my_algs('ecoclin'), [], []);
res.coas = err_res.mean_err.unshuffled;


%% now shuffle post

tr_dir_bins_intrial = tr_dir_bins .* (-1).^within_trial;

if with_conv
    C = tau_conv(S, tau, fps);
else
    C = S;
end
C_shuf = fshuffle(C, tr_dir_bins_intrial, false);

data_tensor = DecodeTensor.construct_tensor(C_shuf, tr_bins, opt.n_bins, tr_s, tr_e);
d.data_tensor = data_tensor;
d.tr_dir = tr_dir;
d.opt = opt;

err_res = DecodeTensor.decode_all(d.data_tensor, d.tr_dir, bin_width, my_algs('ecoclin'), [], []);
res.soc = err_res.mean_err.unshuffled;

tr_dir_bins_intrial = tr_dir_bins .* (-1).^within_trial;
if with_conv
    C = tau_conv(S, tau, fps);
else
    C = S;
end
C_shuf = fshuffle(C, tr_dir_bins_intrial, true);
data_tensor = DecodeTensor.construct_tensor(C_shuf, tr_bins, opt.n_bins, tr_s, tr_e);
d.data_tensor = data_tensor;
d.tr_dir = tr_dir;
d.opt = opt;

err_res = DecodeTensor.decode_all(d.data_tensor, d.tr_dir, bin_width, my_algs('ecoclin'), [], []);
res.asoc = err_res.mean_err.unshuffled;

end

function C = tau_conv(S, ~, ~)
T = 1.5;
Fs = 20;
gamma_sample_times = linspace(0,7.5,T * Fs + 1);
gamma_sample_times = gamma_sample_times(2:end);
gamma = @(t) t.*exp(1-t);
transient_profile = gamma(gamma_sample_times);
transient_profile = transient_profile(:);
C = conv2(S, transient_profile, 'full');
C = C(1:size(S,1),:);
end