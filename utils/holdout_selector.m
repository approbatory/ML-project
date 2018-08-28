function varargout = holdout_selector(holdout_fraction, varargin)
all_cells = all(cellfun(@iscell, varargin));
all_nums = all(cellfun(@isnumeric, varargin));
assert(all_cells || all_nums,...
    'variables to divide must all be cells (indicating trials) or numeric arrays');

assert(isscalar(holdout_fraction) && holdout_fraction > 0 && holdout_fraction < 1, 'holdout_fraction must be a number between 0 and 1');

assert(std(cellfun(@(x) size(x,1), varargin))==0, 'all variables must have the same length on the first axis');
num_samples = size(varargin{1},1);

test_subset = randperm(num_samples)./num_samples <= holdout_fraction;
train_subset = ~test_subset;

varargout = cell(1,2*numel(varargin));
if all_cells
    post_process = @cell2mat;
else
    post_process = @(x)x;
end
for i = 1:numel(varargin)
    varargout{2*i-1} = post_process(varargin{i}(train_subset));
    varargout{2*i} = post_process(varargin{i}(test_subset));
end