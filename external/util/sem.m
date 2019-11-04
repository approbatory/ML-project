function s = sem(x)
% s = sem(x)
%
% Computes the standard error of the mean of the data in the vector x.  
if (isvector(x))
  s = std(x)/sqrt(numel(x));
else
  s = std(x,[],1)./sqrt(size(x,1));
end