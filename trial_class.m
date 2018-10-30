%% PLS of neurons, regressing on fw/bw mean bin position
q_ = o.bin_data(false, false);
mean_bin_X = q_.bin_X.mean;
fw_bindist = mean_bin_X(1:2:end,:) ./ sum(mean_bin_X(1:2:end,:),1);
bw_bindist = mean_bin_X(2:2:end,:) ./ sum(mean_bin_X(2:2:end,:),1);
fw_mb = (1:o.opt.n_bins)*fw_bindist;
bw_mb = (1:o.opt.n_bins)*bw_bindist;
[~, ~, score_,~,~,~,~,stats] = plsregress(o.data.X.fast.', [fw_mb.', bw_mb.'], 2);
neuron_origin = -mean(o.data.X.fast.') * stats.W;
figure; 
subplot(1,2,1); scatter(score_(:,1), score_(:,2), 4, fw_mb);
hold on; scatter(neuron_origin(1), neuron_origin(2), 20, 'r');
subplot(1,2,2); scatter(score_(:,1), score_(:,2), 4, bw_mb);
hold on; scatter(neuron_origin(1), neuron_origin(2), 20, 'r');
suptitle('using PLS on [fw\_mb,bw\_mb] : characterizing neurons');

[~, score_] = pca(o.data.X.fast.', 'NumComponents', 2);
figure; 
subplot(1,2,1); scatter(score_(:,1), score_(:,2), 4, fw_mb);
subplot(1,2,2); scatter(score_(:,1), score_(:,2), 4, bw_mb);
suptitle('using pca : characterizing neurons');
%%
msk = [0; o.data.mask.fast];
trial_start = find(diff(msk) == 1);
trial_end = find(diff(msk) == -1);

i_convert = cumsum(o.data.mask.fast);
num_trials = numel(trial_start);

%% gather stats about each trial
%then display pca of the stats vector, try to detect abnormal trials
trial_feature = zeros(num_trials, o.opt.n_bins, o.data.total_neurons);
trial_dir = zeros(1, num_trials);
for t_i = 1:num_trials
    sel = i_convert(trial_start(t_i):trial_end(t_i));
    sel_save{t_i} = sel;
    trial_X = o.data.X.fast(sel,:);
    trial_ks = o.data.y.ks(sel);
    trial_ks(mod(trial_ks,2)==1) = trial_ks(mod(trial_ks,2)==1) + 1;
    trial_ks = trial_ks/2;
    for b = 1:o.opt.n_bins
        trial_feature(t_i, b, :) = mean(trial_X(trial_ks==b,:));
    end
    trial_dir(t_i) = mean(o.data.y.direction(sel));
end
%%
trial_feature(isnan(trial_feature)) = 0;
trial_feature_flat = reshape(trial_feature, [], o.opt.n_bins*o.data.total_neurons);
[coeff_, score_, latent_] = pca(trial_feature_flat, 'NumComponents', 2);
[coeff_fw, score_fw, latent_fw] = pca(trial_feature_flat(trial_dir==1,:), 'NumComponents', 2);
[coeff_bw, score_bw, latent_bw] = pca(trial_feature_flat(trial_dir==2,:), 'NumComponents', 2);

figure; scatter(score_(:,1), score_(:,2), 4, trial_end - trial_start);
hold on; trial_origin = -mean(trial_feature_flat)*coeff_;
scatter(trial_origin(1), trial_origin(2), 20, 'r', 'filled');
title('Trial-level PCA, color is length');

trial_fw_bw_dim = score_(:,1);
trial_coloring_nonsplit = score_(:,2);
trial_coloring = zeros(1,num_trials);
trial_coloring(trial_dir == 1) = score_fw(:,1);
trial_coloring(trial_dir == 2) = score_bw(:,1);

for t_i = 1:num_trials
    coloring(sel_save{t_i}) = trial_coloring(t_i);
    coloring_nonsplit(sel_save{t_i}) = trial_coloring_nonsplit(t_i);
    coloring_fwbw(sel_save{t_i}) = trial_fw_bw_dim(t_i);
end
%% coloring based on bin mean trial PCA

[XL, ~, act,~,~,~,~,stats] = plsregress(o.data.X.fast, o.data.y.scaled, 2);
%figure; scatter(act(:,1), act(:,2), 4, coloring);
%%
figure; subplot(2,1,1);
scatter(act(:,1), act(:,2), 4, coloring_fwbw);
subplot(2,1,2);
scatter(act(:,1), act(:,2), 4, coloring_nonsplit);
suptitle('trials characterized with PC1&2');
%% coloring based on trial length
trial_length = trial_end - trial_start;
zlength = abs(zscore(trial_length));
for t_i = 1:num_trials
    coloring(sel_save{t_i}) = trial_length(t_i);
end
figure; scatter(act(:,1), act(:,2), 10, coloring, 'filled');
title('Showing trails in PLS colored by length');
%c_ = load('/home/omer/ML-project/plotting/my_colormap_eigen.mat');
%colormap(c_.c);