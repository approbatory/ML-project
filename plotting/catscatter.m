function catscatter(x, y, cat, cat_colors)

cat_vals = unique(cat);

for i = 1:numel(cat_vals)
    color = cat_colors(cat_vals(i));
    if iscell(color) && numel(color)==1
        color = color{1};
    end
    filt = strcmp(cat, cat_vals{i});
    plot(x(filt), y(filt), 'o', 'Color', color);
    hold on;
end