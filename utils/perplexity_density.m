function d = perplexity_density(p)

e = dist_entropy(p);

d = exp(e) ./ size(p,1);