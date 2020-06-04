function varargout = signal_noise_overlap_geo(o) %4
n = @(x)x./norm(x);
dm = cellfun(n, o.dmus, 'UniformOutput', false);


cos_overlap = zeros(2*o.K, 2*o.K-1);
for i = 1:2*o.K
    for j = 1:2*o.K-1
        if isempty(dm{j})
            cos_overlap(i,j) = nan;
            continue;
        end
        cos_overlap(i,j) = cos(subspace(o.evecs{i}(:,1:o.N), dm{j}));
    end
end

if nargout == 1
    varargout{1} = cos_overlap;
    return;
end

imagesc(cos_overlap);
line([20 20]+0.5, ylim, 'Color', 'w');
line(xlim, [20 20]+0.5, 'Color', 'w');
xlabel 'Spatial bin of signal direction';
ylabel 'Spatial bin of noise subspace';
colorbar;
title(sprintf('Cos overlap between %dD noise subspaces\nand signal directions', o.N))
axis image;
end