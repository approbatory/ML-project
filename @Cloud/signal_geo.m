function varargout = signal_geo(o) %2
n = @(x)x./norm(x);
dm = cellfun(n, o.dmus, 'UniformOutput', false);

cos_overlap = zeros(2*o.K-1);
for i = 1:2*o.K-1
    if isempty(dm{i}), continue; end
    for j = 1:2*o.K-1
        if isempty(dm{j}), continue; end
        cos_overlap(i,j) = dm{i}'*dm{j};
    end
end

if nargout == 1
    varargout{1} = cos_overlap;
    return;
end

imagesc(cos_overlap);
xlabel 'Spatial bin'
ylabel 'Spatial bin'
colormap(gca, bluewhitered);
colorbar;
title 'Cos overlap between signal directions'
axis image
end