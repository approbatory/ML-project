function e = dist_entropy(p)

p = p ./ sum(p);

e = p .* log(p);
e(p==0) = 0;

e = -sum(e);