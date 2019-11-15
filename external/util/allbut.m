function Y = allbut(X,i)
% Y = allbut(X,i)
%
% Given a vector X and an integer I, Y = X with the I'th element
% removed. If I<1 or I>numel(X), an out of range error is issued.

if (numel(i)~=1)
  error('Expected I to be a scalar.');
elseif (i<1 | i>numel(X))
  error('Specified index is out of range.');
end

Y = X([1:i-1 i+1:end]);
