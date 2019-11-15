% N_samp = 10000;
% at_zero = randn(N_samp,1);
% at_one = randn(N_samp,1)+2;
%
% X_data = [at_zero; at_one];
% y_data = [zeros(N_samp,1); ones(N_samp,1)];
%
% my_alg = my_algs('lda');
% my_model = my_alg.train(X_data, y_data);
%%

%loading data
save_file = '../linear_track/pablo_data_ds.mat';

%load the data into a ds array
loaded_res = load(save_file);
my_ds = loaded_res.my_ds;

%%
mean_dp2 = @(mdl) sum(sum(mdl.mahal(mdl.Mu)))./(size(mdl.Mu,1).^2-size(mdl.Mu,1));

%arbitrarily choose ix = 18
ix = 18;
ds = my_ds(ix);
midLen = 120 * 10/12;
num_bins = 20;
num_folds = 10;

my_X = ds.trials.traces.';
my_y = ds.trials.centroids;
[sel_fw, sel_bw] = select_directions(my_y); %not using old_method
%just take forwards
my_X = my_X(sel_fw,:);
my_y = my_y(sel_fw,:);
%%
my_binner = @(y) gen_place_bins(y, num_bins, midLen);


my_alg = my_algs('lda');

bin_to_compare = 2;
my_ks = my_binner(my_y);
two_bin_filter = (my_ks == 1) | (my_ks == bin_to_compare);
my_X = my_X(two_bin_filter, :); %added to cut
my_ks = my_ks(two_bin_filter); %added to cut
shuf_X = shuffle(my_X, my_ks);
my_model = my_alg.train(my_X, my_ks);
full_dp2 = mean_dp2(my_model);

part_dp2 = @(n_cells, my_X) mean_dp2(my_alg.train(my_X(:, randperm(size(my_X,2)) <= n_cells), my_ks));
part_single_dp2 = @(n_cells, my_X, my_ks) bw_1n2(my_alg.train(my_X(:, randperm(size(my_X,2)) <= n_cells), my_ks));
%%
num_reps = 20;
delta_neu = 10;
tot_cells = ds.num_cells;
cell_nums = 1:delta_neu:tot_cells;
if cell_nums(end) ~= tot_cells
    cell_nums = [cell_nums tot_cells];
end

dp2_per_cell = zeros(1,numel(cell_nums));
dp2_per_cell_shuf = zeros(1,numel(cell_nums));
for c_ix = numel(cell_nums):-1:1
    num_neu = cell_nums(c_ix);
    for reps = 1:num_reps
        dp2_per_cell(reps,c_ix) = part_single_dp2(num_neu, my_X, my_ks);%part_dp2(i, my_X);
        dp2_per_cell_shuf(reps,c_ix) = part_single_dp2(num_neu, shuf_X, my_ks);%part_dp2(i, shuf_X);
    end
    fprintf('%d ', c_ix);
    if mod(c_ix,30) == 0
        fprintf('\n');
    end
end

fprintf('\n');


%%
%figure; plot(dp2_per_cell); hold on; plot(dp2_per_cell_shuf); legend orig shuf
figure;
hold on;
errb = @(x) std(x)./sqrt(size(x,1));
errorbar(cell_nums, mean(dp2_per_cell), errb(dp2_per_cell));
errorbar(cell_nums, mean(dp2_per_cell_shuf), errb(dp2_per_cell_shuf));

fitline = plot(fishinfo_fit(cell_nums, mean(dp2_per_cell)));
hold on;
fitline_shuf = plot(fishinfo_fit(cell_nums, mean(dp2_per_cell_shuf)));

xlabel('Number of cells');
%ylabel(sprintf('(d'')^2 between bin 1 and %d', bin_to_compare));
ylabel('(d'')^2');
title(sprintf('(d'')^2 between bin 1 and %d vs. Number of cells used', bin_to_compare));

plt = Plot();
fitline.Color = plt.Colors{1};
fitline_shuf.Color = plt.Colors{2};
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.Legend = {'unshuffled', 'shuffled'};
plt.LegendLoc = 'Best';
plt.ShowBox = 'off';
plt.FontSize = 18;
plt.LineStyle = {'none','none','-'};
%%
dp_by_dist = zeros(num_reps,num_bins -1);
dp_by_dist_shuf = zeros(num_reps,num_bins - 1);
for b_ix = 2:num_bins
    for reps = 1:num_reps
    my_X = ds.trials.traces.';
    my_y = ds.trials.centroids;
    [sel_fw, sel_bw] = select_directions(my_y, false); %not using old_method
    %just take forwards
    my_X = my_X(sel_fw,:);
    my_y = my_y(sel_fw,:);
    
    my_binner = @(y) gen_place_bins(y, num_bins, midLen);
    
    
    my_alg = my_algs('lda');
    
    
    my_ks = my_binner(my_y);
    two_bin_filter = (my_ks == 1) | (my_ks == b_ix);
    my_X = my_X(two_bin_filter, :); %added to cut
    my_ks = my_ks(two_bin_filter); %added to cut
    shuf_X = shuffle(my_X, my_ks);
    
    rand_filt = randperm(size(my_X,1))./size(my_X,1) < 0.9;
    
    dp_by_dist(reps, b_ix-1) = part_single_dp2(ds.num_cells, my_X(rand_filt,:), my_ks(rand_filt));
    dp_by_dist_shuf(reps, b_ix-1) = part_single_dp2(ds.num_cells, shuf_X(rand_filt,:), my_ks(rand_filt));
    
    end
    fprintf('%d ', b_ix);
end
fprintf('\n');
%%
figure;
%plot(1:num_bins-1, dp_by_dist./dp_by_dist_shuf, '-o');
hold on;
errorbar(1:num_bins-1, mean(dp_by_dist), errb(dp_by_dist));
errorbar(1:num_bins-1, mean(dp_by_dist_shuf), errb(dp_by_dist_shuf));
legend unshuffled shuffled
title('(d'')^2 between bin 1 and bins 2 through 20');
xlabel('bin distance (bins)');
ylabel('(d'')^2');

dp_rat = dp_by_dist ./ dp_by_dist_shuf;
figure;
errorbar(1:num_bins-1, mean(dp_rat), errb(dp_rat));
title('Ratio of (d'')^2 / (d''_{shuffle})^2');
xlabel('bin distance (bins)');
ylabel('(d'')^2 as fraction of (d''_{shuffle})^2');
%%
function dp2 = bw_1n2(mdl)
res = mdl.mahal(mdl.Mu(1,:));
dp2 = res(2);
end