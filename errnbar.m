function errnbar(y,e)
if isvector(y) || isvector(e)
    y = y(:);
    e = e(:);
end
bar(y);
hold on;
[d,a] = size(y);
ds = 1:d;
as = 1:a;
x = ds' + 1/(a+1.5)*(as - mean(as));
errorbar(x, y, e, '.k');
end