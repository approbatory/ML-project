function output = daysets_end_decoding(daysets, varargin)
p = inputParser;
p.addRequired('daysets', @iscell);
p.addParameter('matched_cells', [], @iscell);
p.addParameter('alg', my_algs('linsvm', 0.1), @isstruct);
p.addParameter('points', 0:0.1:1, @(x) isnumeric(x) && isvector(x) && (min(x) >= 0) && (max(x) <= 1));
p.addParameter('filling', 'copy_zeroed', @ischar);
p.addParameter('par_loops', 64, @(x) isscalar(x) && isnumeric(x));
p.parse(daysets, varargin{:});


points = p.Results.points;
alg = p.Results.alg;
PAR_LOOPS = p.Results.par_loops;
matched_cells = p.Results.matched_cells;
%loop over mice
res = cell(1,numel(daysets));
for m_ix =  1:numel(daysets)
    %loop over days
    for d_ix = 1:numel(daysets{m_ix})
        %loop over points
        if isempty(daysets{m_ix}(d_ix).changing)
            continue;
        end
        fprintf('On m_ix%d d_ix%d\n', m_ix, d_ix);
        %daysets{m_ix}(d_ix).res.points = points;
        res{m_ix}(d_ix).points = points;
        ds = quick_ds(fullfile(daysets{m_ix}(d_ix).directory, daysets{m_ix}(d_ix).day), 'deprobe', 'nocells');
        for j = 1:numel(points)
            point = points(j);
            [X, ks] = ds_dataset(ds, 'selection', point,...
                'filling', p.Results.filling, 'trials', strcmp({ds.trials.start}, daysets{m_ix}(d_ix).changing),...
                'target', {ds.trials.end});
            if nargin == 2
                if size(matched_cells{m_ix},2) ~= numel(daysets{m_ix})%~isfield(daysets{m_ix}(d_ix).res, 'matched_cells')
                    error('Malformed matched cells. %s', daysets{m_ix}(d_ix).day);
                end
                %X = X(:,daysets{m_ix}(d_ix).res.matched_cells); %only use matched cells
                X = X(:, matched_cells{m_ix}(:, d_ix));
            end
            res{m_ix}(d_ix).baseline_err = min(mean(strcmp(ks, 'north')), mean(strcmp(ks, 'south')));
            [res{m_ix}(d_ix).train_err(j,:), res{m_ix}(d_ix).test_err(j,:)] = evaluate_alg(alg,...
                X, strcmp(ks, 'north'), 'par_loops', PAR_LOOPS);
        end
    end
end
output.was_matched = ~isempty(matched_cells);
output.res = res;
end