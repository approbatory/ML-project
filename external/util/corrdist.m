function d = corrdist(X,Y)
% d = corrdist(X,Y)
%
% Computes the correlation distances between the columns of X and Y
% and returns them in the vector d.

if (size(X) ~= size(Y))
  error('X and Y must be the same size.');
end

X = bsxfun(@minus, X, mean(X,1));
Y = bsxfun(@minus, Y, mean(Y,1));

d = 1 - sum(X.*Y,1)./sqrt(sum(X.*X,1))./sqrt(sum(Y.*Y,1));