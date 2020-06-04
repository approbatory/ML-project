function [res, res_sem] = per_sess(o, varname)

res_list = o.mouse_all_sess(varname, []);

notnan = ~isnan(res_list);

res = squeeze(nanmean(res_list, 2));
res_sem = squeeze(nanstd(res_list, 0, 2) ./ sqrt(sum(notnan,2)));