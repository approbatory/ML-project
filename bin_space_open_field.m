function [bins, D] = bin_space_open_field(X, Y)
%assumes position is already preprocessed, rotated, rescaled
if nargin == 1
    Y = X(:,2); X = X(:,1);
end

X = (X - min(X))/(max(X) - min(X));
Y = (Y - min(Y))/(max(Y) - min(Y));
X(X==1) = 1-eps;
Y(Y==1) = 1-eps;
Nx = 7; Ny = 9;
Bx = floor(X * Nx);
By = floor(Y * Ny);
bins = 1 + Bx + Nx*By;

xC = @(B) mod(B - 1, Nx);
yC = @(B) floor((B - 1)/Nx);
dist = @(x,y,x0,y0) sqrt((x-x0)^2 + (y-y0)^2);
D = zeros(Nx*Ny);
for i = 1:Nx*Ny
    for j = 1:Nx*Ny
        x = xC(i); y = yC(i);
        x0 = xC(j); y0 = yC(j);
        D(i,j) = dist(x,y,x0,y0);
    end
end
end