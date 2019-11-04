function Y = mapu(f, X, varargin)
% Y = mapu(f,X, [opts for f])
%
% Applies the function f to each of the columns of X and returns the
% output in the matrix Y. Assumes that the outputs are of uniform
% size.
%
% Y(:,k) = f(X(:,k),varargin{:});
%
% Optional arguments are passed to f.

Y = f(X(:,1),varargin{:});
Y = repmat(Y,1,size(X,2));
for i = 2:size(X,2)
  Y(:,i) = f(X(:,i),varargin{:});
end


    
