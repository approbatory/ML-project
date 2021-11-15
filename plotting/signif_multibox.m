function signif_multibox(X_cell, labels, pairs, top_vals, gap, varargin)
multibox(X_cell, labels, varargin{:});

tonum = containers.Map;
for i = 1:numel(labels)
    tonum(labels{i}) = i;
end

for i = 1:numel(pairs)
    pairs{i} = cellfun(@(x)tonum(x), pairs{i});
    p1 = pairs{i}(1);
    p2 = pairs{i}(2);
    %top_val = max(quantile(X_cell{p1}, 0.75), quantile(X_cell{p2}, 0.75)) + 1;
    top_val = top_vals(i);
    line([p1 p2], [top_val top_val], 'Color', 'k');
    text(mean([p1 p2]), top_val + gap,...
        Utils.pstar(ranksum(X_cell{p1}, X_cell{p2})),...
        'HorizontalAlignment', 'center');
end