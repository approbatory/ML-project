function res = decoder_by_bins_compute(idx, n_reps)

if ~exist('n_reps', 'var')
    n_reps = 20;
end

sm = SessManager;
d = sm.cons_usable(idx);

res.mouse_name = d.mouse_name;
res.source_path = d.source_path;
res.opt = d.opt;

fprintf('Computing RMS error with different sized ensembles\n');

d_neu = 30;
n_neu = size(d.data_tensor, 1);

n_size = unique([1, d_neu:d_neu:n_neu, n_neu]);
num_n_size = numel(n_size);

[me_bins, me_bins_shuf] = deal(zeros(num_n_size, n_reps, d.opt.n_bins));
t_ = tic;
WaitMessage = parfor_wait(n_reps*num_n_size);
for k = 1:num_n_size
    for j = 1:n_reps
    %parfor j = 1:n_reps
        [~, ~, ps, ks, ~] = d.basic_decode(false,n_size(k),[]);
        [~, ~, ps_shuf, ks_shuf, ~] = d.basic_decode(true,n_size(k),[]);
        toc(t_)

        me_bins(k,j,:) = binwise_error(ks, ps, d.opt);
        me_bins_shuf(k,j,:) = binwise_error(ks_shuf, ps_shuf, d.opt);
        WaitMessage.Send;
    end
    fprintf('Finished %d size (%d of %d)\n', n_size(k), k, num_n_size);
end


res.n_size = n_size;
res.me_bins = me_bins;
res.me_bins_shuf = me_bins_shuf;

res.me_bins_max = squeeze(me_bins(end,:,:));
res.me_bins_shuf_max = squeeze(me_bins_shuf(end,:,:));
end

function me_bins = binwise_error(ks, ps, opt)
ks_s = ceil(ks/2);
ps_s = ceil(ps/2);

me_bins = zeros(1, opt.n_bins);
for i = 1:opt.n_bins
    k = ks_s(ks_s == i);
    p = ps_s(ks_s == i);
    
    %me_bins(i) = mean(abs(k - p)) * opt.bin_width;
    me_bins(i) = sqrt(mean((k - p).^2)) * opt.bin_width;
end
end