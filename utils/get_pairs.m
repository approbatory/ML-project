function p = get_pairs(X)
[n,m] = size(X);
assert(n == m, 'not square');

remove_inds = sub2ind([n,m], 1:n, 1:m);
p = X(:);
p(remove_inds) = [];
end