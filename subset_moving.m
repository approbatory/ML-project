function eval_over = subset_moving(ds)
if ds.num_trials == 1
    eval_over(1:ds.full_num_frames) = true;
    return;
end

reindexed_trial_indices = zeros(size(ds.trial_indices));
closing = 0;
for i = 1:size(ds.trial_indices,1)
    r = ds.trial_indices(i,:);
    reindexed_trial_indices(i,:) = r - r(1) + closing + 1;
    closing = reindexed_trial_indices(i,end);
end

eval_over = false(reindexed_trial_indices(end,end),1);
for r = reindexed_trial_indices'
    eval_over(r(2):r(3)) = true;
end
end