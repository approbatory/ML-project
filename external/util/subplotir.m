function inds = subplotir(nr,nc,r,c)
% inds = subplotir(nr,nc,r,c)
%
% Returns the indices for the subplot covering the rows in R and the
% columns in C, in the NR x NC array of subplots.

inds = bsxfun(@plus, nc*(r(:)-1),c(:)');
inds = inds(:);