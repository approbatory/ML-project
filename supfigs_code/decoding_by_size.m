function res_out = decoding_by_size(usable_index, ntype)
if ~exist('ntype', 'var')
    ntype = 'events_transient';
end

sm = SessManager;
d = sm.cons_usable(usable_index, true);
d = DecodeTensor(d, ntype);

alg = my_algs('ecoclin');

n_cells = size(d.data_tensor,1);

d_neu = 30;
n_neu = unique([(d_neu:d_neu:n_cells) n_cells]);
NN = numel(n_neu);

n_reps = 80;

[mse,mse_s,mse_h,mse_hs] = deal(zeros(n_reps, NN));

WaitMessage = parfor_wait(n_reps*NN);
parfor r_i = 1:n_reps
    for i = 1:NN
        
        res = DecodeTensor.decode_all(d.data_tensor, d.tr_dir, d.opt.bin_width, alg, n_neu(i), [], false);
        res_h = DecodeTensor.decode_all(d.data_tensor, d.tr_dir, d.opt.bin_width, alg, n_neu(i), [], true);
        
        mse(r_i, i) = res.MSE.unshuffled;
        mse_s(r_i, i) = res.MSE.shuffled;
        
        mse_h(r_i, i) = res_h.MSE.unshuffled;
        mse_hs(r_i, i) = res_h.MSE.shuffled;
        
        WaitMessage.Send;
    end
end

res_out.mse = mse;
res_out.mse_s = mse_s;
res_out.mse_h = mse_h;
res_out.mse_hs = mse_hs;
res_out.n_neu = n_neu;
res_out.source_path = d.source_path;
res_out.mouse_name = d.mouse_name;
res_out.usable_index = usable_index;
res_out.ntype = ntype;