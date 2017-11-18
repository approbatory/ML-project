function [ds_deprobed] = de_probe(ds)
%DE_PROBE Summary of this function goes here
%   Detailed explanation goes here
n_trials = length(ds.trials);
mask = false(1,n_trials);
for i = 1:n_trials
    mask(i) = strcmp(ds.trials(i).end, 'north') || strcmp(ds.trials(i).end, 'south');
end
ds_deprobed = ds;
ds_deprobed.trials = ds_deprobed.trials(mask);
ds_deprobed.trial_indices = ds_deprobed.trial_indices(mask,:);
ds_deprobed.num_trials = 100;
end

