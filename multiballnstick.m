function [p_val, h] = multiballnstick(labels, ys, errs, varargin)
p = inputParser;
p.addOptional('coloring', {}, @iscell);
p.addOptional('capsize', 1, @isscalar);
p.addOptional('jitter', 0.3, @isscalar);
p.addParameter('medline_color', 'b', @ischar);
p.addParameter('medline_color_alt', 'r', @ischar);
p.addParameter('logscale', false, @islogical);
p.parse(varargin{:});
%shape of all inputs should be 2 x n
n = size(labels, 2);

j_ = p.Results.jitter;%0.15; %jitter magnitude
c_ = p.Results.capsize; %capsize
for i = 1:n
    n_points = numel(ys{1,i});
    %jitter = j_*(rand(1,n_points)-0.5);
    jitter = j_*linspace(-0.5, 0.5, n_points);
    if isempty(p.Results.coloring)
        errorbar((i-1) + repmat([0.85;1.15],1,n_points) + repmat(jitter, 2, 1), [ys{1,i}(:)';ys{2,i}(:)'],...
            [errs{1,i}(:)';errs{2,i}(:)'], '-k', 'Capsize', c_);
        hold on;
        
        errorbar(i+jitter-0.15, ys{1,i}(:)', errs{1,i}(:)', '.b', 'Capsize', c_);
        
        errorbar(0.3+i+jitter-0.15, ys{2,i}(:)', errs{2,i}(:)', '.r', 'Capsize', c_);
    else
        assert(numel(p.Results.coloring) == n, 'coloring cell size must match data size');
        errorbar((i-1) + repmat([0.85;1.15],1,n_points) + repmat(jitter, 2, 1), [ys{1,i}(:)';ys{2,i}(:)'],...
            [errs{1,i}(:)';errs{2,i}(:)'], '-', 'Capsize', c_, 'Color', p.Results.coloring{i});
        hold on;
    end
end

y_pooled_1 = [ys{1,:}];
y_pooled_2 = [ys{2,:}];

l_ = refline(0, median(y_pooled_1)); l_.Color = p.Results.medline_color; l_.LineStyle = ':';
l_ = refline(0, median(y_pooled_2)); l_.Color = p.Results.medline_color_alt; l_.LineStyle = ':';

[p_val, h] = signrank(y_pooled_1, y_pooled_2);
signif_text = Utils.pstar(p_val);

if p.Results.logscale
    set(gca, 'YScale', 'log');
end

yl_ = ylim;
xl_ = xlim;
if p.Results.logscale
    text(mean(xl_), 10.^(log10(yl_(1)) + 0.8*diff(log10(yl_))), signif_text, 'HorizontalAlignment', 'center');
else
    text(mean(xl_), yl_(1) + 0.8*diff(yl_), signif_text, 'HorizontalAlignment', 'center');
end

set(gca, 'XTick', 1:n);
set(gca, 'XTickLabel', labels(1,:));
xlim([0,n+1]);
end