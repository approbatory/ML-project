%function [p, h] = ballnstick(lab1, lab2, y1, y2, e1, e2)
function [p, h] = ballnstick(lab1, lab2, y1, y2, varargin)
p = inputParser;
p.addOptional('e1', []);
p.addOptional('e2', []);
p.addOptional('scatter', false, @islogical);
p.addOptional('coloring', {}, @iscell);
p.addOptional('err_scale', 0.5, @isnumeric);
p.parse(varargin{:});
e1 = p.Results.e1;
e2 = p.Results.e2;

if ~isempty(e1) && ~isempty(e2) && (length(e1) == length(e2))
    n = numel(y1);
    assert(all([numel(y1) numel(y2) numel(e1) numel(e2)] == n), 'sizes must be the same');
    hold on;
    
    if p.Results.scatter
        if isempty(p.Results.coloring)
            scatter([ones(1,n), 2*ones(1,n)]+randn(1,2*n)*0.05,...
                [y1(:)', y2(:)'], p.Results.err_scale./[e1(:)', e2(:)'], 'filled');
        else
            scatter([ones(1,n), 2*ones(1,n)]+randn(1,2*n)*0.05,...
                [y1(:)', y2(:)'], p.Results.err_scale./[e1(:)', e2(:)'],...
                categorical([p.Results.coloring p.Results.coloring]), 'filled');
        end
    else
        errorbar(repmat([1;2], 1, n) + randn(2,n)*0.05, [y1(:)' ; y2(:)'], [e1(:)' ; e2(:)'], '-', 'Capsize', 3);
    end
elseif isempty(e1) && isempty(e2)
    n = numel(y1);
    assert(numel(y1) == numel(y2), 'sizes must be the same');
    plot(repmat([1;2], 1, n) + randn(2,n)*0.05, [y1(:)' ; y2(:)']);
else
    error('e1 and e2 must be errors for y1 and y2');
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