function multibox(X_cell, labels, varargin)

assert(iscell(X_cell));
assert(iscell(labels));
assert(numel(X_cell) == numel(labels));

X_cell = each(@(x)x(:), X_cell);
X_cell = X_cell(:);
num_elems = each(@numel, X_cell);

labels = labels(:);
labels = each(@(l, n) repmat({l}, n, 1), labels, num_elems);
labels = cat(1, labels{:});

boxplot(cell2mat(X_cell), labels, varargin{:});
