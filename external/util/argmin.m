function i = argmin(x,varargin)
% i = argmin(x,varargin)
%
% Finds i such that x(i) = nanmin(x(i)).
%
% An optional argument gives the maximum number of results to return.
if (isempty(varargin))
  i = find(x == nanmin(x));
else
  i = find(x == nanmin(x), varargin{1});
end