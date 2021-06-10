function barerr(X, Y, err)
if isvector(Y) && isvector(err)
    Y = Y(:);
    err = err(:);
end

bar(X, Y);
hold on;
C = bar_center_locations(Y);
errorbar(C, Y, err, 'LineStyle', 'none', 'Color', 'k');
end

function x = bar_center_locations(y)
[d,a] = size(y);
ds = 1:d;
as = 1:a;
x = ds' + 1/(a+1.5)*(as - (a+1)/2);
end
