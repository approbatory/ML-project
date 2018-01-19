function X = ds_dataset(ds, varargin)
p = inputParser;

defaultCombined = true;

defaultSelection = 'all';

defaultFilling = 'copy';
validFilling = {'copy', 'box', 'binary'};
checkFilling = @(x) any(validatestring(x, validFilling));


p.addRequired('ds', @isstruct);
p.addParameter('combined', defaultCombined, @islogical);
p.addParameter('selection', defaultSelection, @checkSelection);
p.addParameter('filling', defaultFilling, checkFilling);

p.parse(ds, varargin{:});

ds = p.Results.ds;
combined = p.Results.combined;
selection = p.Results.selection;
filling = p.Results.filling;

X = gen_place_decoding_X(ds);
if strcmp(filling, 'binary')
    X = cellfun(@(x) x~=0, X, 'UniformOutput', false);
elseif strcmp(filling, 'copy')
    X = cellfun(@(x,y) (x~=0).*(y'), X, {ds.trials.traces}',...
        'UniformOutput', false);
end

if isnumeric(selection)
    [frame_of_interest, ~] = find_frame_at_pos(ds, selection);
    X = cellfun(@(x,i) x(i,:), X, num2cell(frame_of_interest),...
        'UniformOutput', false);
end

if combined
    X = cell2mat(X);
end
end

function res = checkSelection(s)
if ischar(s) && strcmp(s,'all')
    res = true;
    return;
end
if isnumeric(s) && isscalar(s)
    res = true;
    return;
end
res = false;
end