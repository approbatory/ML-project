function [res, res_sem] = all_med_bins(o, varname, restrict)
if ~exist('restrict', 'var')
    restrict = [];
end

[res, res_sem] = o.mouse_med_bins(varname, restrict);