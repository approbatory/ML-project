% to run all:
% S = dir('full_records'); for fs = S', if ~fs.isdir, o = Analyzer.recreate(fullfile(fs.folder, fs.name)); trial_class; end, end; close all; clear; disp('done');
printit = false;
if printit
    mouse_id = split(o.res.source,'/');
    mouse_id = mouse_id{end};
    mouse_id = split(mouse_id, '-');
    mouse_id = [mouse_id{1} mouse_id{2}];
    
    fig_dir = ['graphs2/trial_class/large/' mouse_id];
    if ~exist(fig_dir, 'dir')
        mkdir(fig_dir);
    end
end

%% f_mean_activity_rasters
q_ = o.bin_data(false, false);
mean_bin_X = q_.bin_X.mean;
fw_bindist = mean_bin_X(1:2:end,:) ./ sum(mean_bin_X(1:2:end,:),1);
bw_bindist = mean_bin_X(2:2:end,:) ./ sum(mean_bin_X(2:2:end,:),1);
fw_mb = (1:o.opt.n_bins)*fw_bindist;
bw_mb = (1:o.opt.n_bins)*bw_bindist;
[~, fw_ord] = sort(fw_mb);
[~, bw_ord] = sort(bw_mb);
f_mean_activity_rasters = ...
    figure('Units', 'inches', 'Position',[0 0 10 6]);
subplot(1,2,1); imagesc(mean_bin_X(1:2:end, fw_ord).');
colorbar;
title(sprintf('Mean neural activity\n(forward motion)'));
ylabel('Neuron (ordered by position)');
xlabel('Position bin');
set(gca, 'TickDir', 'out', 'Box', 'off');
subplot(1,2,2); imagesc(mean_bin_X(2:2:end, bw_ord).'); colorbar;
title(sprintf('Mean neural activity\n(backward motion)'));
set(gca, 'TickDir', 'out', 'Box', 'off');
c_ = load('/home/omer/ML-project/plotting/my_colormap_eigen.mat');
colormap(c_.c);
%% f_neuron_pls
[~, ~, score_,~,~,~,~,stats] = plsregress(o.data.X.fast.', [fw_mb.', bw_mb.'], 2);
neuron_origin = -mean(o.data.X.fast.') * stats.W;
f_neuron_pls = figure; 
subplot(1,2,1); scatter(score_(:,1), score_(:,2), 6, fw_mb, 'filled');
hold on; scatter(neuron_origin(1), neuron_origin(2), 20, 'r');
axis equal; xlabel PLS1; ylabel PLS2;
subplot(1,2,2); scatter(score_(:,1), score_(:,2), 6, bw_mb, 'filled');
hold on; scatter(neuron_origin(1), neuron_origin(2), 20, 'r');
axis equal; xlabel PLS1; ylabel PLS2;
suptitle(sprintf('Neuron-level PLS on average bins\nForward and backward directions'));

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
sel_save = cell(1, num_trials);
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
    trial_dir(t_i) = round(mean(o.data.y.direction(sel)));
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

coloring = zeros(size(o.data.X.fast,1),1);
coloring_nonsplit = coloring;
coloring_fwbw = coloring;
for t_i = 1:num_trials
    coloring(sel_save{t_i}) = trial_coloring(t_i);
    coloring_nonsplit(sel_save{t_i}) = trial_coloring_nonsplit(t_i);
    coloring_fwbw(sel_save{t_i}) = trial_fw_bw_dim(t_i);
end
%% f_position_direction_pls

%[XL, ~, act,~,~,~,~,stats] = plsregress(o.data.X.fast, [o.data.y.scaled, o.data.y.direction], 2);
[XL, ~, act,~,~,~,~,stats] = plsregress(o.data.X.fast, zscore([o.data.y.scaled, o.data.y.direction]), 2);
%figure; scatter(act(:,1), act(:,2), 4, coloring);
f_position_direction_pls = figure('Units', 'inches', 'Position', [2.8750 4.4688 11.9792 5.7083]);
subplot(1,3,1);
scatter(act(:,1), act(:,2), 4, o.data.y.direction); axis tight equal;
xlabel PLS1; ylabel PLS2;
title('Colored by direction');
subplot(1,3,2);
scatter(act(:,1), act(:,2), 4, o.data.y.scaled); axis tight equal;
title('Colored by track coord.');
colormap jet
subplot(1,3,3);
scatter(act(:,1), act(:,2), 40, 'black', 'filled', 'MarkerFaceAlpha', 0.015);
origin_ = -mean(o.data.X.fast) * stats.W;
hold on; l_ = scatter(origin_(1), origin_(2), 40, 'red');
axis tight equal;
title('Point density');
legend(l_, 'Origin');
suptitle('Neural activity PLS on direction and track coordinate (z-scored)');
%% f_trialpc_coloring
f_trialpc_coloring = figure('Units', 'inches', 'Position', [3.8646 3.1042 11.3229 5.6354]);
subplot(1,2,1);
scatter(act(:,1), act(:,2), 4, coloring_fwbw); axis tight equal;
xlabel PLS1; ylabel PLS2; title 'Colored by trial-PC1';
subplot(1,2,2);
scatter(act(:,1), act(:,2), 4, coloring_nonsplit); axis tight equal;
xlabel PLS1; ylabel PLS2; title 'Colored by trial-PC2';
suptitle('Trials characterized with trial-PC1&2');
%% f_length_coloring
trial_length = trial_end - trial_start;
zlength = abs(zscore(trial_length));
for t_i = 1:num_trials
    coloring(sel_save{t_i}) = trial_length(t_i);
end
f_length_coloring = figure('Units', 'inches', 'Position', [7.0729 4.2188 5.9479 5.9583]);
scatter(act(:,1), act(:,2), 15, coloring, 'filled', 'MarkerFaceAlpha', 0.2);
axis tight equal;
colorbar;
xlabel PLS1; ylabel PLS2;
title('Colored by trial length (20Hz frames)');
%% f_pls_errors_coloring
total_preds = mean([o.res.unshuf.te_pred{end,:}],2);
decoding_error_coloring = zeros(num_trials,1);
for t_i = 1:num_trials
    decoding_error_coloring(t_i) = mean(abs(o.data.y.ks(sel_save{t_i}) - total_preds(sel_save{t_i})))/2;
end
pointwise_errors = abs(o.data.y.ks - total_preds)/2;

f_pls_errors_coloring = figure('Units', 'inches', 'Position', [3.0208 4.0625 11.4896 5.4583]);
subplot(1,2,1);
scatter(act(pointwise_errors<=1,1), act(pointwise_errors<=1,2),...
    4*(1+pointwise_errors(pointwise_errors<=1)), pointwise_errors(pointwise_errors<=1), 'filled');
axis tight equal;
xl_ = xlim;
yl_ = ylim;
colorbar;
title(sprintf('Bin prediction errors <= 1\n(%.2fcm per bin)',diff(o.data.y.centers(2:3))));
xlabel PLS1; ylabel PLS2;
subplot(1,2,2);
scatter(act(pointwise_errors>1,1), act(pointwise_errors>1,2),...
    4*pointwise_errors(pointwise_errors>1), pointwise_errors(pointwise_errors>1), 'filled');
axis equal;
colormap jet;
colorbar;
xlim(xl_); ylim(yl_);
title('Bin prediction errors > 1');
xlabel PLS1; ylabel PLS2;

%% f_pc2_error
f_pc2_error = figure;
scatter(score_(:,2), decoding_error_coloring);
xlabel 'Trial PC2'; ylabel(sprintf('Mean decoding error per trial\n(%.2fcm per bin)', diff(o.data.y.centers(2:3))));
title('Trial PC2 and mean decoding error');
%% f_trial_pca
f_trial_pca = figure;
scatter(score_(:,1), score_(:,2), 10, trial_dir, 'filled');
hold on; trial_origin = -mean(trial_feature_flat)*coeff_;
l_ = scatter(trial_origin(1), trial_origin(2), 100, 'r');
legend(l_, 'Origin');
title('Trial-level concatenated PCA, colored by direction');
xlabel PC1; ylabel PC2;
colormap redblue
%% f_trial_index_coloring
t_ix_coloring = zeros(1,size(act,1));
for t_i = 1:num_trials
    t_ix_coloring(sel_save{t_i}) = t_i;
end
f_trial_index_coloring = figure('Units', 'inches', 'Position', [0.4062 4.4479 19.3438 4.5833]);
subplot(1,4,4);
scatter(act(:,1), act(:,2), 10, t_ix_coloring, 'filled', 'MarkerFaceAlpha', 0.1);
axis tight equal;
xlabel PLS1; ylabel PLS2;
title('Colored by trial index');
subplot(1,4,1:3);
scatter(1:size(act,1), o.data.y.scaled, 6, t_ix_coloring, 'filled');
xlabel 'Frame (20Hz)'; ylabel 'Position (cm)';
title('Track coordinate');
colormap hsv

%% f_speed_coloring
speed = ([0; diff(o.data.y.raw.full)]);
speed = speed(o.data.mask.fast);
f_speed_coloring = figure('Units', 'inches', 'Position', [0.4062 4.4479 5.5438 4.5833]);
sc_ = speed ./ prctile(speed, 99);
scatter(act(:,1), act(:,2), 20, [max(0,sc_),0*sc_,max(0,-sc_)], 'filled', 'MarkerFaceAlpha', 0.1);
axis tight equal;
colormap redblue
%set(gca, 'Color', 'k');
title('Colored by velocity');
xlabel PLS1; ylabel PLS2;
%% TCA
tr_ft = permute(trial_feature, [3 2 1]);
%n_dim_ = 100;
%err_cp = zeros(1, n_dim_);
%err_cp(1) = 1;
%for dim_ = 2:n_dim_
%    model = cp_nnals(tensor(tr_ft), dim_);
%    temp_ = tr_ft - full(model);
%    err_cp(dim_) = norm(temp_(:))/norm(tr_ft(:));
%    fprintf('Done dim_ = %d\terr_cp = %f\n', dim_, err_cp(dim_));
%end
model = cp_nnals(tensor(tr_ft), 2);
figure; visualize_neuron_ktensor(model, struct('start', trial_dir), 'start', 'sort once');
figure; subplot(2,1,1);
scatter(model.U{1}(:,1), model.U{1}(:,2), 4, fw_mb);
subplot(2,1,2);
scatter(model.U{1}(:,1), model.U{1}(:,2), 4, bw_mb);
figure; scatter(model.U{2}(:,1), model.U{2}(:,2), 20, 1:20);
colormap jet;
figure; scatter(model.U{3}(:,1), model.U{3}(:,2), 20, -trial_dir, 'filled');
colormap redblue

%% printing section
if printit
    printer = @(name) print('-dpng', fullfile(fig_dir, name));
    
    figure(f_mean_activity_rasters);
    printer('mean_activity_rasters.png');
    
    figure(f_neuron_pls);
    printer('neuron_pls.png');
    
    figure(f_position_direction_pls);
    printer('position_direction_pls.png');
    
    figure(f_trialpc_coloring);
    printer('trialpc_coloring.png');
    
    figure(f_length_coloring);
    printer('length_coloring.png');
    
    figure(f_pls_errors_coloring);
    printer('pls_errors_coloring.png');
    
    figure(f_pc2_error);
    printer('pc2_error.png');
    
    figure(f_trial_pca);
    printer('trial_pca.png');
    
    figure(f_trial_index_coloring);
    printer('trial_index_coloring.png');
    
    figure(f_speed_coloring);
    printer('speed_coloring.png');
end