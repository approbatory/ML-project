function [X, Y] = plot_distr(ys, color, x_offset)

num_groups = length(ys);

% Compute the desired percentiles of each group
quartiles = cellfun(@(y) prctile(y,[10 25 50 75 90]), ys, 'UniformOutput', false);
quartiles = cell2mat(quartiles);

% Plot the 25-75 range
x = 1:num_groups;
X = kron(x, [1 1 NaN]);
Y = [quartiles(:,[2 4]) NaN(num_groups,1)]';
plot(X + x_offset, Y(:), 'LineWidth', 1, 'Color', color);
hold on;

% % Plot the 10-90 range
% Y = [quartiles(:,[1 5]) NaN(num_groups,1)]';
% plot(X + x_offset, Y(:), 'LineWidth', 1, 'Color', color);

% Plot the median
plot(x + x_offset, quartiles(:,3), '.', 'MarkerSize', 16, 'Color', color);
