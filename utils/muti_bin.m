function MI = muti_bin(X,Y, K, counts_y)
N = length(X);
if (length(X) ~= length(Y)) || ~isvector(X) || ~isvector(Y)
    error('must be vectors of equal length');
end

X = X(:); Y = Y(:);

if nargin == 2
    K = max(Y);
    TAB = tabulate(Y);
    counts_y = TAB(:,2);
end

counts_one = sparse(Y(X==1),1,1,K,1);
%counts_one = zeros(K,1);
%TAB = tabulate(Y(X==1));
%TAB2 = TAB(:,2);
%counts_one(1:length(TAB2)) = TAB2;
counts_zero = counts_y - counts_one;
counts = [counts_zero, counts_one]';

Pxy = counts / N;
Px = sum(Pxy,2); Py = sum(Pxy,1);
%MI = Pxy .* (log(Pxy) - log(Px) - log(Py));
MI = Pxy .* (log(Pxy./(Px.*Py)));
MI(isinf(MI) | isnan(MI)) = 0;
MI = sum(sum(MI));
MI = MI/log(2);

end