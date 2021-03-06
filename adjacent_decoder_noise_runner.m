function res = adjacent_decoder_noise_runner(index, no_decoding, n_type)
if ~exist('no_decoding', 'var')
    no_decoding = false;
end
if ~exist('n_type', 'var')
    n_type = 'rawTraces';
end
tic
%d_ = DecodeTensor.cons_filt(index);
sm = SessManager;
o_ = sm.cons_usable(index, true);
d_ = DecodeTensor(o_, n_type);
n_reps = 80;
[N, K, T] = size(d_.data_tensor);
n_sizes = unique([2, 10:10:N, N]);
for j = 1:n_reps
    for i = 1:numel(n_sizes)
        %[res.nv{j,i}, res.nv_s{j,i}, res.nv_d{j,i},...
        %    res.m2{j,i}, res.m2_s{j,i}, res.m2_d{j,i}, res.md2{j,i}] = ...
        %    DecodeTensor.adjacent_decoder_noise(d_.data_tensor, d_.tr_dir, n_sizes(i), []);
        res.results_table{j,i} = ...
            DecodeTensor.adjacent_metrics(d_.data_tensor, d_.tr_dir, n_sizes(i), [], no_decoding);
        progressbar(j/n_reps, i/numel(n_sizes));
    end
end
res.source = d_.source_path;
res.filt_index = index;
res.n_sizes = n_sizes;
toc
