function [bins, D] = bin_space(X,Y)
%assumes position is already preprocessed, rotated, rescaled
if nargin == 1
    Y = X(:,2); X = X(:,1);
end
bins = zeros(size(X));
X(X==1) = 1-eps;
Y(Y==1) = 1-eps;
N = 10;
X_bins = mod(floor(N*X),N)+1;
Y_bins = mod(floor(N*Y),N)+1;
tb = (Y>(1-X)) == (Y>X);
bins(tb) = Y_bins(tb);
bins(~tb) = N+X_bins(~tb);

D = dist_matrix(N);
end


function D = dist_matrix(N)
%N even
D = -ones(2*N);
for i = 1:N
    for j = 1:N
        D(i,j) = abs(i-j);
    end
end

for i = N+1:2*N
    for j = N+1:2*N
        D(i,j) = abs(i-j);
    end
end

for i = N+1:2*N
    for j = 1:N
        mid = [N/2, N/2+1];
        d1 = min(abs(i-N - mid));
        d2 = min(abs(j   - mid));
        D(i,j) = d1 + d2 + 1;
        D(j,i) = D(i,j);
    end
end
end