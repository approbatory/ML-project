function MI = muti_bin(iX,Y, K, counts_y)
N = numel(Y);
%if (length(X) ~= length(Y)) || ~isvector(X) || ~isvector(Y)
%    error('must be vectors of equal length');
%end

%X = X(:); 
Y = Y(:);

if nargin == 2
    error('supply all arguments');
    K = max(Y);
    TAB = tabulate(Y);
    counts_y = TAB(:,2);
end

counts_one = sparse(Y(iX),1,1,K,1);
counts_zero = counts_y - counts_one;
counts = [counts_zero, counts_one]';

Pxy = counts / N;
Px = sum(Pxy,2); Py = sum(Pxy,1);
MI = Pxy .* (log(Pxy./(Px.*Py)));
MI(isinf(MI) | isnan(MI)) = 0;
MI = sum(MI(:));
MI = MI/log(2);

end