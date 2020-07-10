function [cell_mean, cell_sem] = cell_mean_sem(C, dim)


N = size(C, dim);

sz = size(C);
out_sz = sz;
out_sz(dim) = 1;

[cell_mean, cell_sem] = cell(out_sz);

out_numel = numel(cell_mean);

for i = 1:out_numel
    <> use vec_ind2sub %%TODO
%% TODO
end



