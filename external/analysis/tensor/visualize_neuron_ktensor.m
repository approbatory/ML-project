function [Ax, G, BigAx] = visualize_neuron_ktensor(decomp, meta, trialcolor, neuron_sort_method)
% VISUALIZE_NEURON_KTENSOR, wraps VISUALIZE_KTENSOR for multi-day data
%
%     [Ax, BigAx, FigHandle] = VISUALIZE_NEURON_KTENSOR(X, meta, trialcolor, neuron_sort_method)
%
%     Parameters
%     ----------
%     trialcolor : 'start', 'error', 'end', etc.

if nargin == 2
    neuron_sort_method = 'sort once';
end


% color for the trials
c = cell(1,3);
if nargin >=3
    c{3} = categorical_colors(meta.(trialcolor));
end

% the neuron order is arbitrary, so permute/sort it along the first factor
decomp.u{1} = decomp.u{1} .* repmat(decomp.lambda', size(decomp.u{1}, 1), 1);
decomp.lambda = ones(size(decomp.lambda));

switch neuron_sort_method
case 'sort once'
    [~, n_idx] = sort(sum(abs(decomp.u{1}), 2), 'descend');
    decomp.u{1} = decomp.u{1}(n_idx, :);
case 'sort each'
    for r = 1:size(decomp.u{1}, 2)
       [~, n_idx] = sort(abs(decomp.u{1}(:,r)), 'descend');
       decomp.u{1}(:,r) = decomp.u{1}(n_idx, r);
    end
otherwise
    error('Sorting option for neurons not recognized.')
end

% bar plot for neuron factors
% line plot for within trial factors
% line and scatter plot for across trial factors
plt = {'bar', 'line', 'scatter'};

% plot titles
nm = {'neuron factors', 'time factors',...
      sprintf('trial factors (%s)', trialcolor)};

% formatting
lspc = {[], '-r', '-k'};
lw = [0 2 1];

[Ax, G, BigAx] = visualize_ktensor(decomp, 'c', c, 'plots', plt, ...
                  'linewidth', lw, 'link_yax', [true, true, true], ...
                  'title', nm, 'linespec', lspc);


set(Ax(1:end-1,:), 'YTickLabel', [])
set(Ax, 'YAxisLocation', 'left')
set(Ax(end, :), 'TickDir', 'out')

w_factor = [0.02, 0.0, 0.02];
[M, N] = size(Ax);
for m = 1:M
    for n = 1:N
        pos = get(Ax(m,n), 'Position');
        pos(1) = pos(1) - w_factor(n);
        pos(3) = pos(3) + 2*w_factor(n);
        set(Ax(m,n), 'Position', pos)
    end
end

for r = 1:size(Ax,1)
  set(get(Ax(r,1), 'ylabel'), 'String', sprintf('r = %g ',r), 'Rotation', 0, 'HorizontalAlignment', 'right');
end

set(get(Ax(end,1), 'xlabel'), 'String', ['method = ' neuron_sort_method]);
set(Ax(end,:), 'YAxisLocation', 'right')
set(get(Ax(end,1), 'ylabel'), 'Position', get(get(Ax(1,1), 'ylabel'), 'Position'));
