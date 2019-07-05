function multiballnstick(labels, ys, errs)
%shape of all inputs should be 2 x n
n = size(labels, 2);

j_ = 0.3;%0.15; %jitter magnitude
c_ = 1; %capsize
for i = 1:n
    n_points = numel(ys{1,i});
    %jitter = j_*(rand(1,n_points)-0.5);
    jitter = j_*linspace(-0.5, 0.5, n_points);
    errorbar((i-1) + repmat([0.85;1.15],1,n_points) + repmat(jitter, 2, 1), [ys{1,i}(:)';ys{2,i}(:)'],...
        [errs{1,i}(:)';errs{2,i}(:)'], '-k', 'Capsize', c_);
    hold on;
    
    errorbar(i+jitter-0.15, ys{1,i}(:)', errs{1,i}(:)', '.b', 'Capsize', c_);
    
    errorbar(0.3+i+jitter-0.15, ys{2,i}(:)', errs{2,i}(:)', '.r', 'Capsize', c_);
    
end

set(gca, 'XTick', 1:n);
set(gca, 'XTickLabel', labels(1,:));
xlim([0,n+1]);
end