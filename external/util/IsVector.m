function b = IsVector(X)
% b = IsVector(X)
%
% Returns TRUE if X is a vector, FALSE otherwise.

b = isnumeric(X) & numel(size(X))==2 & min(size(X))==1;
