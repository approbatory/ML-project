function q = perm_inv(p)

q = 1:numel(p);
q(p) = q;