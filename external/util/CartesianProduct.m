function [P,I] = CartesianProduct(X, varargin)
% [P,I] = CartesianProduct(X)
%
% Given an M x N generalized matrix X, returns a Prod_i
% (size(X(:,i))) by N matrix P whose rows contain the values of the
% cartesian product of the columns of X. The matrix I is of the
% same size as P and returns the indices into the columns of X used
% for each element of the cartesian product.
%
% Altenatively, X can be a vector and an optional argument can
% specify the number of repetitions N of X, so that the cartesian
% product Pi_N X is computed. An additionaly argument set to
% 'Unique' will only return the product corresponding to unique
% sets of indices.

onlyUnique = 0;

if (numel(varargin)>2)
  error('At most 2 optional arguments expected.');
elseif (numel(varargin)>0)
  if (~isvector(X))
    error('When two arguments are provided, the first must be a vector.');
  end
  numRepetitions = varargin{1};
  X = repmat(X(:),1,numRepetitions);
  if (numel(varargin)>1)
    onlyUnique = isequal(upper(varargin{2}), 'UNIQUE');
  end
end

if (~iscell(X))
  X = map(@Identity, X);
end

numRows = cellfun(@numel, X);

I = cell2mat(arrayfun(@(x) moddecomp(x, numRows)', 0:prod(numRows)-1, 'UniformOutput', false))'+1;

if (onlyUnique)
  I = unique(sort(I,2),'rows');
end

P = cell2mat(mapi(@(inds,n) Columnize(X{n}(inds)), I));