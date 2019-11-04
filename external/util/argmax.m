function i = argmax(x,varargin)
% i = argmax(x,...)
%
% Finds i such that x(i) = max(x(i)). An optional argument can
% supply the maximum number of search results. By default all are
% returned.

if (~isempty(varargin))
  i = find(x == max(x), varargin{1});
else
  i = find(x == max(x));
end
