function dist = mean_dist(k,p)
[~,D] = bin_space([],[]);
dist_func = D(sub2ind(size(D), k, p));
dist = mean(dist_func);
end