function Y = mapi(f, X, varargin)
% Y = mapi(f,X, [opts for f])
%
% Indexed map. Applies the function f to each of the columns of X
% and returns the output in a cell array Y, i.e.
%
% Y{k} = f(X(:,k),k,varargin{:});
%
% Note that the index k is passed as the second argument to f.
% Optional arguments are passed to f.

Y = cell(1,size(X,2));

for i = 1:size(X,2)
  Y{i} = f(X(:,i),i,varargin{:});
end


    
