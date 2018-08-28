function varargout = kfold_selector(k, k_i, varargin)
all_cells = all(cellfun(@iscell, varargin));
all_nums = all(cellfun(@isnumeric, varargin));
assert(all_cells || all_nums,...
    'variables to divide must all be cells (indicating trials) or numeric arrays');

assert(isscalar(k) && (round(k) == k) && (k > 0), 'k must be a positive integer');
assert(isscalar(k_i) && (round(k_i) == k_i) && (k_i > 0) && (k_i <= k), 'k_i must be a positive integer in the range 1:k');

%assert(std(cellfun(@(x) size(x,1), varargin))==0, 'all variables must have the same length on the first axis');
%num_samples = size(varargin{1},1);

kfold_divisions = @(num_samples) floor((0:num_samples-1)./num_samples.*k) + 1;
%test_subset = kfold_divisions == k_i;
%train_subset = ~test_subset;

varargout = cell(1,2*numel(varargin));
if all_cells
    post_process = @cell2mat;
else
    post_process = @(x)x;
end
for i = 1:numel(varargin)
    divs = kfold_divisions(size(varargin{i},1));
    test_subset = divs == k_i;
    train_subset = ~test_subset;
    varargout{2*i-1} = post_process(varargin{i}(train_subset,:));
    varargout{2*i} = post_process(varargin{i}(test_subset,:));
end