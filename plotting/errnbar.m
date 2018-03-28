function errnbar(y,e, overlay)
overlay =  exist('overlay', 'var');
if isvector(y) || isvector(e)
    y = y(:);
    e = e(:);
end
if ~overlay
    bar(y);
end
hold on;
[d,a] = size(y);
ds = 1:d;
as = 1:a;
x = ds' + 1/(a+1.5)*(as - mean(as));
errorbar(x, y, e, '.k');
if overlay
    plot(x, y, 'o');
end
end