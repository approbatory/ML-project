function bins = bin_space(X,Y)
%assumes position is already preprocessed, rotated, rescaled
if nargin == 1
    Y = X(:,2); X = X(:,1);
end
X(X==1) = 1-eps;
Y(Y==1) = 1-eps;
N = 10;
X_bins = mod(floor(N*X),N)+1;
Y_bins = mod(floor(N*Y),N)+1;
tb = (Y>(1-X)) == (Y>X);
bins(tb) = Y_bins(tb);
bins(~tb) = N+X_bins(~tb);
end