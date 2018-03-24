function [X, ks, eval_metric] = ds_dataset(ds, varargin)
p = inputParser;

defaultCombined = true;

defaultSelection = 'all';

defaultFilling = 'copy';
validFilling = {'copy', 'box', 'binary', 'traces', 'copy_zeroed'};
checkFilling = @(x) any(validatestring(x, validFilling));

%Trial selection parameters:
defaultTrials = 'all';
checkTrials = @(t) (ischar(t) && strcmp(t,'all')) ||...
    (islogical(t) && isvector(t) && (numel(t) == ds.num_trials));

%Learning target parameter:
defaultTarget = 'position bin';
checkTarget = @(t) (ischar(t) && strcmp(t, 'position bin')) ||...
    (iscell(t) && isvector(t) && (numel(t) == ds.num_trials));

defaultSparsify = true;

defaultOpenField = false;

p.addRequired('ds', @isstruct);
p.addParameter('combined', defaultCombined, @islogical);
p.addParameter('selection', defaultSelection, @checkSelection);
p.addParameter('filling', defaultFilling, checkFilling);
p.addParameter('trials', defaultTrials, checkTrials);
p.addParameter('target', defaultTarget, checkTarget);
p.addParameter('sparsify', defaultSparsify, @islogical);
p.addParameter('openfield', defaultOpenField, @islogical);

p.parse(ds, varargin{:});

ds = p.Results.ds;
combined = p.Results.combined;
selection = p.Results.selection;
filling = p.Results.filling;
trials = p.Results.trials;
target = p.Results.target;
sparsify = p.Results.sparsify;
openfield = p.Results.openfield;

X = gen_place_decoding_X(ds);
if strcmp(filling, 'binary')
    X = cellfun(@(x) double(x~=0), X, 'UniformOutput', false);
elseif strcmp(filling, 'copy')
    X = cellfun(@(x,y) (x~=0).*(y'), X, {ds.trials.traces}',...
        'UniformOutput', false);
elseif strcmp(filling, 'copy_zeroed')
    X = gen_place_decoding_X(ds, true);
elseif strcmp(filling, 'traces')
    X = cellfun(@(x)x', {ds.trials.traces}', 'UniformOutput', false);
end


if ischar(target) && strcmp(target, 'position bin')
    if openfield
        bin_func = @bin_space_open_field;
    else
        bin_func = @bin_space;
    end
    [~,D] = bin_func([],[]);
    dist_func = @(k,p) D(sub2ind(size(D), k(:), p(:))); %accomodate row/col vecs
    mean_dist = @(k,p) mean(dist_func(k,p));
    eval_metric = mean_dist;
    if openfield
        XY = {ds.trials.centroids};
    else
        XY = preprocess_xy(ds);
    end
    ks = cellfun(bin_func, XY, 'UniformOutput', false);
else
    eval_metric = @(k,p) mean(~cellfun(@isequal, k, p));
    ks = target;
end


if isnumeric(selection)
    [frame_of_interest, ~] = find_frame_at_pos(ds, selection);
    X = cellfun(@(x,i) x(i,:), X, num2cell(frame_of_interest),...
        'UniformOutput', false);
    if ischar(target) && strcmp(target, 'position bin')
        ks = cellfun(@(x,i) x(i), ks, num2cell(frame_of_interest),...
            'UniformOutput', false);
    end
elseif strcmp(selection, 'moving')
    starts = ds.trial_indices(:,2) - ds.trial_indices(:,1) + 1;
    ends = ds.trial_indices(:,3) - ds.trial_indices(:,1);
    X = cellfun(@(x,s,e) x(s:e,:), X, num2cell(starts), num2cell(ends),...
        'UniformOutput', false);
    if ischar(target) && strcmp(target, 'position bin')
        ks = cellfun(@(x,s,e) x(s:e), ks, num2cell(starts), num2cell(ends),...
            'UniformOutput', false);
    end
end


if islogical(trials)
    X = X(trials);
    ks = ks(trials);
end

if combined
    X = cell2mat(X);
    if sparsify
        X = sparse(X);
    end
    if ~iscell(target)
        ks = cell2mat(ks);
    end
elseif sparsify
    X = cellfun(@sparse, X, 'UniformOutput', false);
end
end

function res = checkSelection(s)
if ischar(s) && (strcmp(s,'all') || strcmp(s,'moving'))
    res = true;
    return;
end
if isnumeric(s) && isscalar(s)
    res = true;
    return;
end
res = false;
end