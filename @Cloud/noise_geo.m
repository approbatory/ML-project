function varargout = noise_geo(o) %3
cos_overlap = zeros(2*o.K);
for i = 1:2*o.K
    for j = 1:2*o.K
        cos_overlap(i,j) = mean(cos(subspacea(o.evecs{i}(:,1:o.N), o.evecs{j}(:,1:o.N))));
    end
end

if nargout == 1
    varargout{1} = cos_overlap;
    return;
end

imagesc(cos_overlap);
line([20 20]+0.5, ylim, 'Color', 'w');
line(xlim, [20 20]+0.5, 'Color', 'w');
xlabel 'Spatial bin';
ylabel 'Spatial bin';
colorbar;
title(sprintf('Mean canon corr between\n%dD noise subspaces', o.N));
axis image
end