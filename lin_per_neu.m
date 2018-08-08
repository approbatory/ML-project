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

alg = my_algs('ecoclin');
alg_preshuf = my_algs('ecoclin', 'shuf');
num_samples = 64*8;
midLen = 0.8*120;
num_bins = 20;


my_X = ds.trials.traces.';
my_y = ds.trials.centroids;
[sel_fw, sel_bw] = select_directions(my_y);
%just take forwards
my_X = my_X(sel_fw,:);
my_y = my_y(sel_fw,:);

my_binner = @(y) gen_place_bins(y, num_bins, midLen);

delta_neu = 5;%10;
tot_cells = ds.num_cells;
cell_nums = 1:delta_neu:tot_cells;
tr_err = zeros(num_samples, numel(cell_nums));
te_err = zeros(num_samples, numel(cell_nums));
tr_err_preshuf = zeros(num_samples, numel(cell_nums));
te_err_preshuf = zeros(num_samples, numel(cell_nums));
tr_err_shufboth = zeros(num_samples, numel(cell_nums));
te_err_shufboth = zeros(num_samples, numel(cell_nums));
for c_ix = 1:numel(cell_nums)
    my_ticker = tic;
    num_neu = cell_nums(c_ix);
    
    parfor rep_ix = 1:num_samples
        cell_subset_X = my_X(:,randperm(tot_cells)<=num_neu);
        [tr_err(rep_ix,c_ix), te_err(rep_ix,c_ix)] = evala(alg, cell_subset_X, my_y, my_binner,...
            'split', 'nonlocal', 'verbose', false, 'errfunc', 'RMS');
        [tr_err_preshuf(rep_ix,c_ix), te_err_preshuf(rep_ix,c_ix)] = evala(alg_preshuf, cell_subset_X, my_y, my_binner,...
            'split', 'nonlocal', 'verbose', false, 'errfunc', 'RMS');
        [tr_err_shufboth(rep_ix,c_ix), te_err_shufboth(rep_ix,c_ix)] = evala(alg, cell_subset_X, my_y, my_binner,...
            'split', 'nonlocal', 'verbose', false, 'errfunc', 'RMS', 'shufboth', true);
        fprintf('%d ', rep_ix);
    end
    fprintf('\n');
    toc(my_ticker);
    fprintf('For %d cells: %.2f +- %.2f RMS\t ORIGINAL\n', num_neu, mean(te_err(:,c_ix)), std(te_err(:,c_ix))/sqrt(num_samples));
    fprintf('For %d cells: %.2f +- %.2f RMS\t PRESHUF\n', num_neu, mean(te_err_preshuf(:,c_ix)), std(te_err_preshuf(:,c_ix))/sqrt(num_samples));
    fprintf('For %d cells: %.2f +- %.2f RMS\t SHUFBOTH\n', num_neu, mean(te_err_shufboth(:,c_ix)), std(te_err_shufboth(:,c_ix))/sqrt(num_samples));
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
plt.export(sprintf('graphs/place_decoding/linear_track_hpc/info_per_neuron/lin_track_day%d_info_per_neuron_fixed_shufboth.png', ix));

%% saving data
save(sprintf('records/lin_track_day%d_info_per_neuron_larger_fixed_shufboth_before_split.mat', ix), 'tr_err', 'te_err', 'tr_err_preshuf', 'te_err_preshuf', 'tr_err_shufboth', 'te_err_shufboth');

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
suptitle('[I_0N]/(1+N\epsilon) fit parameters');
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

%% messaging
mailme('Job done', sprintf('The job of lin_per_neu.m for day index %d has finished', ix));