function [ks, centers, dims, scXY, numrowcol]  = gen_place_bins(ds, divs, scale, islinear)
if ~exist('islinear', 'var')
    islinear = false;
end
if isnumeric(ds)
    XY = ds;
else
    if ds.num_trials ~= 1
        error('must have one trial, open field');
    end
    XY = ds.trials.centroids;
end
X = XY(:,1);
if size(XY,2) == 2
    Y = XY(:,2);
end

x_range = max(X) - min(X);
if size(XY,2) == 2
    y_range = max(Y) - min(Y);
end


if (size(XY,2) == 2) && (y_range > x_range)
    temp = X;
    X = Y;
    Y = temp;
    
    temp_r = x_range;
    x_range = y_range;
    y_range = temp_r;
end %make X always longer

X = (X - min(X))./x_range;
X(X==1) = 1-eps;

if size(XY,2) == 2
    Y = (Y - min(Y))./x_range;
    y_frac_range = max(Y) - min(Y);
    Y(Y./y_frac_range == 1) = y_frac_range - eps;
end

Nx = divs; 

if size(XY,2) == 2
    if islinear
        Ny = 1;
    else
        Ny = max(1,floor(y_frac_range * Nx));
    end
end

if size(XY,2) == 2
    numrowcol = [Nx, Ny];
else
    numrowcol = Nx;
end

Bx = floor(X * Nx);
if size(XY,2) == 2
    By = floor(Y ./ y_frac_range * Ny);

    if islinear && false
        Direc = (diff(X) > 0) + 1;
        Direc = [2; Direc];
        ks = sub2ind([2, Nx, Ny], Direc, Bx+1, By+1);
    else
        ks = 1 + Bx + Nx*By;
        %ks = sub2ind([Nx, Ny], Bx+1, By+1);
    end
    dims = [1 ./ Nx, y_frac_range ./ Ny].*scale;
else
    ks = 1 + Bx;
    dims = scale ./ Nx;
end

if size(XY,2) == 2
    centers_func = @(i,j) [i-1/2,j-1/2].*dims;

    centers = zeros(Nx*Ny, 2);
    for k = 1:Nx*Ny
        i = mod((k-1), Nx)+1;
        j = floor((k-1)./Nx)+1;
        centers(k,:) = centers_func(i,j);
    end
else
    centers = ((1:Nx)-(1/2))*dims;
end

scX = X*scale;
if size(XY,2) == 2
    scY = Y*scale;
    scXY = [scX,scY];
else
    scXY = scX;
end