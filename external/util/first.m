function f = first(X,varargin)
% f = first(X,varargin)
%
% Returns the first k elements of X. k is by default 1, but can be
% specified with an optional numeric argument. A second optional
% argument set to 'rand' will return the first k points after random
% shuffle.
if (numel(X)==0)
  error('Input vector is empty.');
end
if (numel(varargin)>=2 && isequal(lower(varargin{2}),'rand'))
  X = X(randperm(numel(X)));
end

if (isempty(varargin))
  f = X(1);
else
  f = X(1:varargin{1});
end