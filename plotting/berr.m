function berr(X, Ys_main, Ys_sec, leg)
if length(X)~=size(Ys_main,1) ||...
        ~isequal(size(Ys_main),size(Ys_sec))
    error(['X must have an entry for each row in Y,'...
        ' Ys_main and Ys_sec must have the same size']);
end

sem = @(x) std(x)/sqrt(length(x));
Ys_main_means = cellfun(@mean, Ys_main);
Ys_main_errbs = cellfun(sem, Ys_main);
Ys_sec_means = cellfun(@mean, Ys_sec);
Ys_sec_errbs = cellfun(sem, Ys_sec);

bar(Ys_main_means);
hold on;
C = bar_center_locations(Ys_main_means);
errorbar(C, Ys_main_means, Ys_main_errbs, '.k');
plot(C, Ys_sec_means, 'o');
errorbar(C, Ys_sec_means, Ys_sec_errbs, '.k');
legend(leg{:});
set(gca, 'XTickLabel', X);
set(gca, 'XTick', 1:length(X));

function x = bar_center_locations(y)
[d,a] = size(y);
ds = 1:d;
as = 1:a;
x = ds' + 1/(a+1.5)*(as - (a+1)/2);
end
end