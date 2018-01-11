function MI = muti(X,Y,Vx,Vy)
N = length(X);
if (length(X) ~= length(Y)) || ~isvector(X) || ~isvector(Y)
    error('must be vectors of equal length');
end

X = X(:); Y = Y(:);

if nargin == 2
    Vx = unique(X);
    Vy = unique(Y);
end
Vx = Vx(:); Vy = Vy(:);

Pxy = double(X' == Vx) * double(Y == Vy') / N;
Px = sum(Pxy,2); Py = sum(Pxy,1);
MI = Pxy .* (log(Pxy) - log(Px) - log(Py));
MI(isinf(MI) | isnan(MI)) = 0;
MI = sum(sum(MI));
MI = MI/log(2);
end