function [res, res_sem] = mouse_med_bins(o, varname, mouse_name)

res_list = o.mouse_all_sess(varname, mouse_name);

%3 - sess, 2 - bins, 1 - dims

res_list = nanmedian(res_list, 2);

res = mean(res_list, 3);
res_sem = std(res_list, 0, 3) ./ sqrt(size(res_list,3));