function plotmat(titl, train_errs, test_errs, algs, dayset, err_lim, self_std)
if ~exist('self_std', 'var')
    self_std = false;
end
if self_std
    errb = @(x)std(x);
else
    errb = @(x)std(x)/sqrt(length(x));
end
errnbar(cellfun(@mean, test_errs), cellfun(errb, test_errs));
errnbar(cellfun(@mean, train_errs), cellfun(errb, train_errs), 'overlay');
set(gca, 'XTickLabel', {dayset.day});
set(gca, 'XTick', 1:numel({dayset.day}));
xtickangle(5);
legend(algs.name);
ylim([0 err_lim]);
ylabel('bin distance');
title(titl);
end