function [res, res_sem] = all_by_bins(o, varname, restrict)
if ~exist('restrict', 'var')
    restrict = [];
end

[res, res_sem] = o.mouse_by_bins(varname, restrict);