function daysets = daysets_end_decoding(daysets, matching)
points = 0:0.1:1;
alg = my_algs('linsvm', 0.1); %0.1 L1 regularization
PAR_LOOPS = 64;
%loop over mice
for m_ix =  1:numel(daysets)
    %loop over days
    for d_ix = 1:numel(daysets{m_ix})
        %loop over points
        if isempty(daysets{m_ix}(d_ix).changing)
            continue;
        end
        fprintf('On m_ix%d d_ix%d\n', m_ix, d_ix);
        daysets{m_ix}(d_ix).res.points = points;
        for j = 1:numel(points)
            point = points(j);
            [X, ks] = ds_dataset(daysets{m_ix}(d_ix).ds, 'selection', point,...
                'filling', 'copy_zeroed', 'trials', strcmp({daysets{m_ix}(d_ix).ds.trials.start}, daysets{m_ix}(d_ix).changing),...
                'target', {daysets{m_ix}(d_ix).ds.trials.end});
            if matching
                if ~isfield(daysets{m_ix}(d_ix).res, 'matched_cells')
                    error('Does not have matched cells. %s', daysets{m_ix}(d_ix).day);
                end
                X = X(:,daysets{m_ix}(d_ix).res.matched_cells); %only use matched cells
            end
            daysets{m_ix}(d_ix).res.baseline_err = min(mean(strcmp(ks, 'north')), mean(strcmp(ks, 'south')));
            [daysets{m_ix}(d_ix).res.train_err(j,:), daysets{m_ix}(d_ix).res.test_err(j,:)] = evaluate_alg(alg,...
                X, strcmp(ks, 'north'), 'par_loops', PAR_LOOPS);
        end
    end
end
end