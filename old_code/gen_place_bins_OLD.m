function [ks, centers, dims, scXY, numrowcol]  = gen_place_bins(ds, divs, scale)
if isnumeric(ds)
    XY = ds;
else
    if ds.num_trials ~= 1
        error('must have one trial, open field');
    end
    XY = ds.trials.centroids;
end
X = XY(:,1); Y = XY(:,2);

x_range = max(X) - min(X);
y_range = max(Y) - min(Y);

if y_range > x_range
    temp = X;
    X = Y;
    Y = temp;
    
    temp_r = x_range;
    x_range = y_range;
    y_range = temp_r;
end %make X always longer

X = (X - min(X))./x_range;
X(X==1) = 1-eps;
Y = (Y - min(Y))./x_range;
y_frac_range = max(Y) - min(Y);
Y(Y./y_frac_range == 1) = y_frac_range - eps;
Nx = divs; Ny = max(1,floor(y_frac_range * Nx));
numrowcol = [Nx, Ny];

Bx = floor(X * Nx);
By = floor(Y ./ y_frac_range * Ny);
ks = 1 + Bx + Nx*By;

dims = [1 ./ Nx, y_frac_range ./ Ny].*scale;

centers_func = @(i,j) [i-1/2,j-1/2].*dims;

centers = zeros(Nx*Ny, 2);
for k = 1:Nx*Ny
    i = mod((k-1), Nx)+1;
    j = floor((k-1)./Nx)+1;
    centers(k,:) = centers_func(i,j);
end

scX = X*scale;
scY = Y*scale;
scXY = [scX,scY];