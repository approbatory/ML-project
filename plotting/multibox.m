function multibox(X_cell, labels, varargin)

assert(iscell(X_cell));
assert(iscell(labels));
assert(numel(X_cell) == numel(labels));

X_cell = cellfun(@(x)x(:), X_cell, 'UniformOutput', false);
X_cell = X_cell(:);
num_elems = cellfun(@numel, X_cell, 'UniformOutput', false);

labels = labels(:);
labels = cellfun(@(l, n) repmat({l}, n, 1), labels, num_elems, 'UniformOutput', false);
labels = cat(1, labels{:});

boxplot(cell2mat(X_cell), labels, varargin{:});
