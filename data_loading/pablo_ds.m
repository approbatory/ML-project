function [ds, recfile] = pablo_ds(tracesEvents)
[ds.full_num_frames, ds.num_cells] = size(tracesEvents.rawTraces);
ds.trial_indices = [1 2 3 ds.full_num_frames];
ds.num_trials = 1;
ds.trials.start = 'east';
ds.trials.goal = 'north';
ds.trials.end = 'north';
ds.trials.time = 1;
ds.trials.correct = 1;
ds.trials.turn = 'right';
ds.trials.centroids = tracesEvents.position(1:ds.full_num_frames,:);
ds.trials.traces = tracesEvents.rawTraces.';

recfile.traces = tracesEvents.rawTraces;
recfile.filters = ones(100,100, ds.num_cells);
recfile.info.type = 'cellmax';
recfile.info.num_pairs = ds.num_cells;
end