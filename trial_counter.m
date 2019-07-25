fid = fopen('trial_counter.txt', 'w');
progressbar('Counting trials');
for i = 1:239
    opt = DecodeTensor.default_opt;
    opt.n_bins = 10;
    try
        d_ = DecodeTensor(DecodeTensor.cons(i,true), [], opt);
        fprintf(fid, '%d\n', size(d_.data_tensor,3));
    catch e
        fprintf(fid, '-1\n');
    end
    progressbar(i/239);
end