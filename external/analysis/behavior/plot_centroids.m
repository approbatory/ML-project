function plot_centroids(ds, trial_inds, pos_frames)

if ~iscell(trial_inds) % Allow for multiple groups of trial_inds
    trial_inds = {trial_inds};
end

group_colors = 'gbr';
num_groups = length(trial_inds);

[height, width] = size(ds.behavior_ref_img);
As = zeros(height, width, ds.num_trials); % Preallocate
a_idx = 0;

centroids = cell(num_groups,1);
for j = 1:num_groups
    if islogical(trial_inds{j})
        trial_inds_j = find(trial_inds{j});
    else
        trial_inds_j = trial_inds{j};
    end
    num_trials_in_j = length(trial_inds_j);
    centroids_j = zeros(num_trials_in_j,2);
    for k = 1:num_trials_in_j
        trial_idx = trial_inds_j(k);
        pos_frame = pos_frames(trial_idx);

        centroids_j(k,:) = ds.trials(trial_idx).centroids(pos_frame,:);
        
        a_idx = a_idx + 1;
        As(:,:,a_idx) = ds.get_behavior_trial_frame(trial_idx, pos_frame);
    end
    centroids{j} = centroids_j;
end
As = As(:,:,1:a_idx);

A = mean(As,3);
imagesc(A);
axis image;
colormap gray;
hold on;
for j = 1:num_groups
    color = group_colors(mod(j,length(group_colors))+1);
    plot(centroids{j}(:,1), centroids{j}(:,2), 'o',...
        'MarkerSize', 8,...
        'MarkerEdgeColor', 'w',...
        'MarkerFaceColor', color);
end
hold off;