function n = normc(X, varargin)
% n = normc(X, varargin);
%
% Takes the p-norm along the columns. p defaults to 2, but can be set
% with an optional second argument.

if (isempty(varargin))
  p = 2;
else  
  p = varargin{1};
  if (p<0)
    error('p must be positive.');
  end
end

if (isinf(p))
  n = max(X);
elseif (p==0)
  n = sum(X~=0,1);
else
  n = sum(abs(X).^p,1).^(1/p);
end