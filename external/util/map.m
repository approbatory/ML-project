function Y = map(f, X, varargin)
% Y = map(f,X, [opts for f])
%
% Applies the function f to each of the columns of X and returns
% the output in a cell array Y, i.e.
%
% Y{k} = f(X(:,k),varargin{:});
%
% Optional arguments are passed to f.

Y = cell(1,size(X,2));

for i = 1:size(X,2)
  Y{i} = f(X(:,i),varargin{:});
end