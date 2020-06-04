function res = sess_med_bins(o, varname, sess_idx)

res = o.sess_by_bins(varname, sess_idx);

res = nanmedian(res, 2);