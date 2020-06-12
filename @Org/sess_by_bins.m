function res = sess_by_bins(o, varname, sess_idx)

V = o.fetch(varname);

V = V{sess_idx};

assert(iscell(V), 'not a cell, maybe not showing all bins');

num_dims = cellfun(@numel, V);
num_dims(num_dims == 0) = [];
n = min(num_dims);

V = cellfun(@(x)sel_(x,n), V, 'UniformOutput', false);
V = cellfun(@(x)x(:), V, 'UniformOutput', false);
arr = cell2mat(V);
%imagesc(arr);

res = arr;


end

function r = sel_(x, n)

if isempty(x)
    r = nan(n,1);
else
    r = x(1:n);
end
end