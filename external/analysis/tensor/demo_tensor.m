% NOTE: Need to pre-load data as a MultiDay object called `md`
[X, meta, neuron_map, trial_map] = export(md);

% TODO: consider standardize before warp?
% X = standardize(X,trial_map);

% preprocessing (timewarp and then standardize)
X = timewarp(X);
% for c = 1:size(X,1)
%     x = X(c,:,:);
%     X(c,:,:) = x ./ (max([1, max(abs(x(:)))]));
% end

% make a scree plot
[cpd_list,rsq] = fit_cpd(X,'min_rank',1);
visualize_rank(cpd_list);

% pick best cpd to analyze further
[~,i] = max(rsq);
cpd = cpd_list(i);
Xest = full(cpd.decomp);
Xest = Xest.data;

% plot single-figure summary of all factors
visualize_neuron_ktensor(cpd.decomp,md,trial_map);

% plot fit across neurons
visualize_fit(X,Xest,1,md,trial_map);

% plot fit across trials
visualize_fit(X,Xest,3,md,trial_map);

% TO SAVE FIGURES TO A DIRECTORY:
% output_dir = '/path/to/destination/'
% visualize_fit(X,Xest,1,md,trial_map,output_dir);
% visualize_fit(X,Xest,3,md,trial_map,output_dir);

% plot residuals
visualize_resids(X,Xest,md,trial_map);
