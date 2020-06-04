function [res, res_sem] = mouse_by_bins(o, varname, mouse_name)

res_list = o.mouse_all_sess(varname, mouse_name);

res = mean(res_list, 3);
res_sem = std(res_list, 0, 3) ./ sqrt(size(res_list,3));