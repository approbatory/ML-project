function [p, h] = ballnstick(lab1, lab2, y1, y2, e1, e2)
if nargin == 6
    n = numel(y1);
    assert(all([numel(y1) numel(y2) numel(e1) numel(e2)] == n), 'sizes must be the same');
    
    errorbar(repmat([1;2], 1, n) + randn(2,n)*0.05, [y1(:)' ; y2(:)'], [e1(:)' ; e2(:)'], '-', 'Capsize', 3);
elseif nargin == 4
    n = numel(y1);
    assert(numel(y1) == numel(y2), 'sizes must be the same');
    plot(repmat([1;2], 1, n) + randn(2,n)*0.05, [y1(:)' ; y2(:)']);
else
    error('nargin must be 4 or 6');
end
set(gca, 'XTick', [1 2]);
set(gca, 'XTickLabel', {lab1, lab2});
xlim([0.5 2.5]);

[p, h] = signrank(y1, y2);

if p < 0.001
    signif_text = '***';
elseif p < 0.01
    signif_text = '**';
elseif p < 0.05
    signif_text = '*';
else
    signif_text = 'n.s';
end
yl_ = ylim;
text(1.5, yl_(1) + 0.8*diff(yl_), signif_text, 'HorizontalAlignment', 'center');
end