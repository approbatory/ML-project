function r = value_rank(x)

r = zeros(size(x));

x_flat = x(:);
[~, ord] = sort(x_flat);
r_flat = perm_inv(ord);

r(1:numel(r)) = r_flat;