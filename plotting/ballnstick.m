%function [p, h] = ballnstick(lab1, lab2, y1, y2, e1, e2)
function [p, h] = ballnstick(lab1, lab2, y1, y2, varargin)
p = inputParser;
p.addOptional('e1', []);
p.addOptional('e2', []);
p.addOptional('scatter', false, @islogical);
p.addOptional('coloring', {}, @iscell);
p.addOptional('jitter', 0.3, @isnumeric);
p.addOptional('err_scale', 0.5, @isnumeric);
p.addOptional('relative_err', false, @islogical);
p.parse(varargin{:});
e1 = p.Results.e1;
e2 = p.Results.e2;
assert(numel(p.Results.coloring)==numel(y1)||isempty(p.Results.coloring), 'coloring cell (numel %d) does not match data size (%d)', numel(p.Results.coloring), numel(y1));

if ~isempty(e1) && ~isempty(e2) && (length(e1) == length(e2))
    n = numel(y1);
    assert(all([numel(y1) numel(y2) numel(e1) numel(e2)] == n), 'sizes must be the same');
    hold on;
    
    if p.Results.scatter
        if p.Results.relative_err
            err_marker = p.Results.err_scale./[e1(:)'./y1(:)', e2(:)'./y2(:)'];
        else
            err_marker = p.Results.err_scale./[e1(:)', e2(:)'];
        end
        if isempty(p.Results.coloring)
            scatter([ones(1,n), 2*ones(1,n)]+randn(1,2*n)*p.Results.jitter,...
                [y1(:)', y2(:)'], err_marker, 'filled');
        else
            scatter([ones(1,n), 2*ones(1,n)]+randn(1,2*n)*p.Results.jitter,...
                [y1(:)', y2(:)'], err_marker,...
                categorical([p.Results.coloring p.Results.coloring]), 'filled');
        end
    else
        if isempty(p.Results.coloring)
            errorbar(repmat([1;2], 1, n) + linspace(-0.5, 0.5, n)*p.Results.jitter, [y1(:)' ; y2(:)'], [e1(:)' ; e2(:)'], '-', 'Capsize', 3);
        else
            offset = linspace(-0.5, 0.5, n);
            for i = 1:n
                hold on;
                errorbar([1;2] + offset(i)*p.Results.jitter, [y1(i) ; y2(i)],...
                    [e1(i) ; e2(i)], '-', 'Capsize', 3, 'Color', p.Results.coloring{i});
            end
        end
    end
elseif isempty(e1) && isempty(e2)
    n = numel(y1);
    assert(numel(y1) == numel(y2), 'sizes must be the same');
    if isempty(p.Results.coloring)
        plot(repmat([1;2], 1, n) + linspace(-0.5, 0.5, n)*p.Results.jitter, [y1(:)' ; y2(:)']);
    else
        offset = linspace(-0.5, 0.5, n);
        for i = 1:n
            hold on;
            plot([1;2] + offset(i)*p.Results.jitter, [y1(i) ; y2(i)],...
                'Color', p.Results.coloring{i});
        end
    end
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