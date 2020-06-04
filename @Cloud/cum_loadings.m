function N_50 = cum_loadings(o, suppress) %1
if ~exist('suppress', 'var')
    suppress = false;
end


%colors = hsv(o.K);
for i = 1:2*o.K-1
    if isempty(o.loadings{i})
        continue;
    end
    cs = cumsum(o.loadings{i}.^2);
    cs_s = cumsum(o.loadings_shuf{i}.^2);
    reorder = @(x) x(randperm(numel(x)));
    cs_ro = cumsum(reorder(o.loadings{i}.^2));
    %plot(cs, 'Color', colors(mod(i-1,o.K)+1,:));
    if ~suppress
        hold on
        plot(cs, 'b');
        plot(cs_s, 'r');
        plot(cs_ro, 'g');
    end
    covers_half = find(cs > 0.5, 1);
    if ~isempty(covers_half)
        ind(i) = covers_half;
    else
        ind(i) = numel(cs);
        fprintf('Not enough components to reach 50%% signal variance (only up to %f). Using maximal number of %d\n', cs(end), ind(i));
    end
end
ind = ind(ind ~= 0);
%figure;
%histogram(ind);
N_50 = round(median(ind));
fprintf('Median index to 50%% is %d\n', N_50);

if ~suppress
    xlabel 'Number of PCs'
    ylabel 'Fraction of \Delta\mu within the subspace'
    text(160, 0.2, 'Real', 'Color', 'b');
    text(160, 0.15, 'Shuffled', 'Color', 'r');
    text(160, 0.1, 'Unordered', 'Color', 'g');
    title 'Showing all spatial bins'
end
end