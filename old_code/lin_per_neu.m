%% lin_per_neu.m


%loading data
save_file = '../linear_track/pablo_data_ds.mat';

%load the data into a ds array
loaded_res = load(save_file);
my_ds = loaded_res.my_ds;

%%
%arbitrarily choose ix = 18
ix = 18;
ds = my_ds(ix);

report_mahal = true;
if report_mahal
    alg = my_algs('lda');
    alg_preshuf = my_algs('dlda');
else
    alg = my_algs('ecoclin');
    alg_preshuf = my_algs('ecoclin', 'shuf');
end

num_samples = 2;%20;%16;%64
midLen = 120 * 10/12;
num_bins = 20;
num_folds = 10;

my_X = ds.trials.traces.';
my_y = ds.trials.centroids;
[sel_fw, sel_bw] = select_directions(my_y, false); %not using old_method
%just take forwards
my_X = my_X(sel_fw,:);
my_y = my_y(sel_fw,:);

my_binner = @(y) gen_place_bins(y, num_bins, midLen);

delta_neu = 10;
tot_cells = ds.num_cells;
cell_nums = 1:delta_neu:tot_cells;
if cell_nums(end) ~= tot_cells
    cell_nums = [cell_nums tot_cells];
end
tr_err = cell(num_samples, numel(cell_nums));
te_err = cell(num_samples, numel(cell_nums));
tr_err_preshuf = cell(num_samples, numel(cell_nums));
te_err_preshuf = cell(num_samples, numel(cell_nums));
tr_err_shufboth = cell(num_samples, numel(cell_nums));
te_err_shufboth = cell(num_samples, numel(cell_nums));
%for c_ix = 1:numel(cell_nums)
for c_ix = numel(cell_nums):-1:1
    my_ticker = tic;
    num_neu = cell_nums(c_ix);
    
    %for rep_ix = 1:num_samples
    parfor rep_ix = 1:num_samples
        cell_subset_X = my_X(:,randperm(tot_cells)<=num_neu);
        [tr_err{rep_ix,c_ix}, te_err{rep_ix,c_ix}, models_ord] = evala(alg, cell_subset_X, my_y, my_binner,...
            'split', 'nonlocal', 'verbose', false, 'errfunc', 'mean_dist', 'kfold', num_folds, 'error_type', 'binwise', 'train_frac', 0.9);
        [tr_err_preshuf{rep_ix,c_ix}, te_err_preshuf{rep_ix,c_ix}, models_preshuf] = evala(alg_preshuf, cell_subset_X, my_y, my_binner,...
            'split', 'nonlocal', 'verbose', false, 'errfunc', 'mean_dist', 'kfold', num_folds, 'error_type', 'binwise', 'train_frac', 0.9);
        [tr_err_shufboth{rep_ix,c_ix}, te_err_shufboth{rep_ix,c_ix}, models_shufboth] = evala(alg, cell_subset_X, my_y, my_binner,...
            'split', 'nonlocal', 'verbose', false, 'errfunc', 'mean_dist', 'shufboth', true, 'kfold', num_folds, 'error_type', 'binwise', 'train_frac', 0.9);
        if report_mahal
            mahal_ord(rep_ix, c_ix) = mean(mean_mahal_dists_calculator(models_ord));
            mahal_preshuf(rep_ix, c_ix) = mean(mean_mahal_dists_calculator(models_preshuf));
            mahal_shufboth(rep_ix, c_ix) = mean(mean_mahal_dists_calculator(models_shufboth));
        end
        
        
        fprintf('%d ', rep_ix);
    end
    
    fprintf('\n');
    toc(my_ticker);
    fprintf('For %d cells: %.2f +- %.2f mean err\t ORIGINAL\n', num_neu, mean(cell2mat(te_err(:,c_ix))), std(cell2mat(te_err(:,c_ix)))/sqrt(num_samples));
    fprintf('For %d cells: %.2f +- %.2f mean err\t PRESHUF\n', num_neu, mean(cell2mat(te_err_preshuf(:,c_ix))), std(cell2mat(te_err_preshuf(:,c_ix)))/sqrt(num_samples));
    fprintf('For %d cells: %.2f +- %.2f mean err\t SHUFBOTH\n', num_neu, mean(cell2mat(te_err_shufboth(:,c_ix))), std(cell2mat(te_err_shufboth(:,c_ix)))/sqrt(num_samples));
    if report_mahal
        fprintf('For %d cells: %.2f +- %.2f mean mahal\t ORIGINAL\n', num_neu, mean(mahal_ord(:,c_ix)), std(mahal_ord(:,c_ix))/sqrt(num_samples));
        fprintf('For %d cells: %.2f +- %.2f mean mahal\t PRESHUF\n', num_neu, mean(mahal_preshuf(:,c_ix)), std(mahal_preshuf(:,c_ix))/sqrt(num_samples));
        fprintf('For %d cells: %.2f +- %.2f mean mahal\t SHUFBOTH\n', num_neu, mean(mahal_shufboth(:,c_ix)), std(mahal_shufboth(:,c_ix))/sqrt(num_samples));
    end
end
tr_err = cell2mat(tr_err); te_err = cell2mat(te_err);
tr_err_preshuf = cell2mat(tr_err_preshuf); te_err_preshuf = cell2mat(te_err_preshuf);
tr_err_shufboth = cell2mat(tr_err_shufboth); te_err_shufboth = cell2mat(te_err_shufboth);


%% saving data
if report_mahal
    save(sprintf('records/lin_track_day_mahal_%d_%s.mat', ix, timestring), 'tr_err', 'te_err', 'tr_err_preshuf', 'te_err_preshuf', 'tr_err_shufboth', 'te_err_shufboth');
else 
    save(sprintf('records/lin_track_day%d_%s.mat', ix, timestring), 'tr_err', 'te_err', 'tr_err_preshuf', 'te_err_preshuf', 'tr_err_shufboth', 'te_err_shufboth');
end
%% messaging
%mailme('Re: Job done', sprintf('The job of lin_per_neu.m for day index %d has finished (not old_method)', ix));
%return;
%% if return_mahal
if report_mahal
    mahal_ord_mean = mean(mahal_ord); mahal_ord_errb = std(mahal_ord)./sqrt(num_samples);
    mahal_preshuf_mean = mean(mahal_preshuf); mahal_preshuf_errb = std(mahal_preshuf)./sqrt(num_samples);
    mahal_shufboth_mean = mean(mahal_shufboth); mahal_shufboth_errb = std(mahal_shufboth)./sqrt(num_samples);
    
    figure;
    hold on;
    errorbar(cell_nums, mahal_ord_mean, mahal_ord_errb, 'b');
    errorbar(cell_nums, mahal_preshuf_mean, mahal_preshuf_errb, 'r');
    errorbar(cell_nums, mahal_shufboth_mean, mahal_shufboth_errb, 'g');
    
    xlabel('Number of cells');
    ylabel('Mean (d'')^2');
    title('Mean (d'')^2 vs. Number of cells used');
    
    plt = Plot();
    plt.XMinorTick = 'off';
    plt.YMinorTick = 'off';
    plt.Legend = {'ordinary', 'trained on shuffle', 'trained & tested on shuffle'};
    plt.LegendLoc = 'Best';
    plt.ShowBox = 'off';
    plt.FontSize = 18;
    return;
end

%% Transform to Fisher information

te_fish = 1./(te_err.^2);
te_fish_mean = mean(te_fish);
te_fish_errb = std(te_fish)./sqrt(num_samples);

te_fish_preshuf = 1./(te_err_preshuf.^2);
te_fish_mean_preshuf = mean(te_fish_preshuf);
te_fish_errb_preshuf = std(te_fish_preshuf)./sqrt(num_samples);

te_fish_shufboth = 1./(te_err_shufboth.^2);
te_fish_mean_shufboth = mean(te_fish_shufboth);
te_fish_errb_shufboth = std(te_fish_shufboth)./sqrt(num_samples);

%%
figure;
hold on;
errorbar(cell_nums, te_fish_mean, te_fish_errb, 'b');
errorbar(cell_nums, te_fish_mean_preshuf, te_fish_errb_preshuf, 'r');
errorbar(cell_nums, te_fish_mean_shufboth, te_fish_errb_shufboth, 'g');

fitline = plot(fishinfo_fit(cell_nums, te_fish_mean));
fitline_preshuf = plot(fishinfo_fit(cell_nums, te_fish_mean_preshuf));
fitline_shufboth = plot(fishinfo_fit(cell_nums, te_fish_mean_shufboth));

xlabel('Number of cells');
ylabel('Fisher information (cm^{-2})');
title('Information vs. Number of cells used');

plt = Plot();
fitline.Color = plt.Colors{1};
fitline_preshuf.Color = plt.Colors{2};
fitline_shufboth.Color = plt.Colors{3};

plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.Legend = {'ordinary', 'trained on shuffle', 'trained & tested on shuffle'};
plt.LegendLoc = 'Best';
plt.ShowBox = 'off';
plt.FontSize = 18;
plt.LineStyle = {'none','none','none', '-'};
%plt.Colors = [0 0 1; 1 0 0; 0 1 0; 0 0 1];
%plt.export(sprintf('graphs/place_decoding/linear_track_hpc/info_per_neuron/lin_track_day%d_info_per_neuron_fixed_shufboth.png', ix));
%% showing as normal mean error
te_err_mean = mean(te_err);
te_err_errb = std(te_err)./sqrt(num_samples);
te_err_mean_preshuf = mean(te_err_preshuf);
te_err_errb_preshuf = std(te_err_preshuf)./sqrt(num_samples);
te_err_mean_shufboth = mean(te_err_shufboth);
te_err_errb_shufboth = std(te_err_shufboth)./sqrt(num_samples);

figure;
hold on;
errorbar(cell_nums, te_err_mean, te_err_errb, 'b');
errorbar(cell_nums, te_err_mean_preshuf, te_err_errb_preshuf, 'r');
errorbar(cell_nums, te_err_mean_shufboth, te_err_errb_shufboth, 'g');

%fitline = plot(fishinfo_fit(cell_nums, te_fish_mean));
%fitline_preshuf = plot(fishinfo_fit(cell_nums, te_fish_mean_preshuf));
%fitline_shufboth = plot(fishinfo_fit(cell_nums, te_fish_mean_shufboth));

xlabel('Number of cells');
ylabel('Mean error (cm)');
title('Error vs. Number of cells used');

plt = Plot();
%fitline.Color = plt.Colors{1};
%fitline_preshuf.Color = plt.Colors{2};
%fitline_shufboth.Color = plt.Colors{3};

plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.Legend = {'ordinary', 'trained on shuffle', 'trained & tested on shuffle'};
plt.LegendLoc = 'Best';
plt.ShowBox = 'off';
plt.FontSize = 18;
%plt.LineStyle = {'none','none','none', '-'};
%plt.Colors = [0 0 1; 1 0 0; 0 1 0; 0 0 1];
%plt.export(sprintf('graphs/place_decoding/linear_track_hpc/info_per_neuron/lin_track_day%d_info_per_neuron_fixed_shufboth.png', ix));


%% fit params shown
ordinary_fit = fishinfo_fit(cell_nums, te_fish_mean);
preshuf_fit = fishinfo_fit(cell_nums, te_fish_mean_preshuf);
shufboth_fit = fishinfo_fit(cell_nums, te_fish_mean_shufboth);
%%

Y_errnbar = [coeffvalues(ordinary_fit);coeffvalues(preshuf_fit);coeffvalues(shufboth_fit)];
E_errnbar = [diff(confint(ordinary_fit));diff(confint(preshuf_fit));diff(confint(shufboth_fit))];
figure;
subplot(2,1,1);
errnbar(Y_errnbar(:,1), E_errnbar(:,1));
set(gca, 'XTickLabels', {'Ordinary','Train on \newline shuffle','Train & test \newline on shuffle'});
ylabel('I_0 (cm^{-2}/neuron)');
%plt_fit1 = Plot();
subplot(2,1,2);
errnbar(Y_errnbar(:,2), E_errnbar(:,2));
set(gca, 'XTickLabels', {'Ordinary','Train on \newline shuffle','Train & test \newline on shuffle'});
ylabel('\epsilon (1/neuron)');
suplabel('[I_0N]/(1+N\epsilon) fit parameters');
%plt_fit2 = Plot();

%%

I_shuffle_mean = te_fish_mean - te_fish_mean_shufboth;
I_shuffle_errb = sqrt(te_fish_errb.^2 + te_fish_errb_shufboth.^2);

I_diag_mean = te_fish_mean - te_fish_mean_preshuf;
I_diag_errb = sqrt(te_fish_errb.^2 + te_fish_errb_preshuf.^2);

figure;
hold on;

errorbar(cell_nums, I_diag_mean, I_diag_errb);
errorbar(cell_nums, I_shuffle_mean, I_shuffle_errb);

xlabel('Number of cells');
ylabel('\Delta Fisher information (cm^{-2})');
title('Information gain compared to shuffled data');
legend '\DeltaI_{diag}' '\DeltaI_{shuffle}' 'Location' 'best'

plt_diff = Plot();
plt_diff.Colors{2} = plt.Colors{3};
plt_diff.Colors{1} = plt.Colors{2};
plt_diff.XMinorTick = 'off';
plt_diff.YMinorTick = 'off';
plt_diff.ShowBox = 'off';

