function P = CartesianProduct1(v, n)
% P = CartesianProduct1(v, n)
%
% Given the vector v, returns a matrix whose rows contain all n-tuples
% of elements of v.

X = cell(n,1);

cmdStr = '[';
for i = 1:n-1
  cmdStr = sprintf('%sX{%d}, ', cmdStr, i);
end
cmdStr = sprintf('%sX{%d}] = ndgrid(v);',cmdStr, n);
eval(cmdStr);

P = zeros(numel(X{1}), n);
for i = 1:n
  P(:,i) = X{i}(:);
end
  

