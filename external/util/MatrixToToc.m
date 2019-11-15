function Y = MatrixToToc(X,d)
% Y = MatrixToToc(X,d)
%
% Converts the matrix X into a toc matrix Y with the
% dimensions d. X is assumed to be in 
% (Trials x Odors x Cells x Values) format.

% Rearrange so that the values for a given trial are along the first dimension.
X = permute(X,[4 1 2 3]);
Y = reshape(X,[],prod(d));
