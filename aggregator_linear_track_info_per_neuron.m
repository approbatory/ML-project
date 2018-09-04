%record_pattern = 'records_sherlock_315/Mouse-2024-20150315_090450-linear-track-TracesAndEvents.mat_decoding_*.mat';
record_pattern = 'records/lin_track_2024_0317_*.mat';
using = 'last'; %or 'all'

%search for all files matching the pattern
S_dir = dir(record_pattern);

switch using
    case 'last'
        [~, order] = sortrows({S_dir.name}.');
        S_dir = S_dir(order(end));
    case 'all'
    otherwise
        error('specify ''last'' or ''all'' for using');
end

total_res = [];
total_res_shuf = [];
for i = 1:numel(S_dir)
    loaded = load(fullfile(S_dir(i).folder, S_dir(i).name));
    total_res = [total_res loaded.res];
    total_res_shuf = [total_res_shuf loaded.res_shuf];
end
%%
[me_m, me_e, n_cells] = fetch_from_rec(total_res, 'mean_err');
[me_sh_m, me_sh_e, n_cells_sh] = fetch_from_rec(total_res_shuf, 'mean_err');

[mse_m, mse_e, ~] = fetch_from_rec(total_res, 'MSE');
[mse_sh_m, mse_sh_e, ~] = fetch_from_rec(total_res_shuf, 'MSE');

[fi_m, fi_e, ~] = fetch_from_rec(total_res, 'fisher_info');
[fi_sh_m, fi_sh_e, ~] = fetch_from_rec(total_res_shuf, 'fisher_info');

%% plot
figure;
errorbar(n_cells, me_m, me_e);
hold on;
errorbar(n_cells_sh, me_sh_m, me_sh_e);
legend unshuffled shuffled
xlabel('Number of cells');
ylabel('Mean error (cm)');
title('Mean decoding error vs. Number of cells used');

figure;
errorbar(n_cells, mse_m, mse_e);
hold on;
errorbar(n_cells_sh, mse_sh_m, mse_sh_e);
legend unshuffled shuffled
xlabel('Number of cells');
ylabel('Mean squared error (cm^2)');
title('Mean squared decoding error vs. Number of cells used');

figure;
errorbar(n_cells, fi_m, fi_e);
hold on;
errorbar(n_cells_sh, fi_sh_m, fi_sh_e);
legend unshuffled shuffled
xlabel('Number of cells');
ylabel('Fisher information (cm^{-2})');
title('Fisher information vs. Number of cells used');

%% fancy fisher information plot
figure;
errorbar(n_cells, fi_m, fi_e);
hold on;
errorbar(n_cells_sh, fi_sh_m, fi_sh_e);

fitline = plot(fishinfo_fit(n_cells, fi_m));
fitline_sh = plot(fishinfo_fit(n_cells_sh, fi_sh_m));

legend unshuffled shuffled Location best
xlabel('Number of cells');
ylabel('Fisher information (cm^{-2})');
title('Fisher information vs. Number of cells used');

plt = Plot();
fitline.Color = plt.Colors{1};
fitline_sh.Color = plt.Colors{2};

plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
%plt.Legend = {'ordinary', 'trained on shuffle', 'trained & tested on shuffle'};
plt.LegendLoc = 'Best';
plt.ShowBox = 'off';
plt.FontSize = 18;
plt.LineStyle = {'none','none', '-'};

%% fit params shown
unshuffled_fit = fishinfo_fit(n_cells, fi_m);
shuffled_fit = fishinfo_fit(n_cells_sh, fi_sh_m);


Y_errnbar = [coeffvalues(unshuffled_fit);coeffvalues(shuffled_fit)];
E_errnbar = [diff(confint(unshuffled_fit));diff(confint(shuffled_fit))];
figure;
subplot(2,1,1);
errnbar(Y_errnbar(:,1), E_errnbar(:,1));
set(gca, 'XTickLabels', {'Unshuffled','Shuffled'});
ylabel('I_0 (cm^{-2}/neuron)');
%plt_fit1 = Plot();
subplot(2,1,2);
errnbar(Y_errnbar(:,2), E_errnbar(:,2));
set(gca, 'XTickLabels', {'Unshuffled','Shuffled'});
ylabel('\epsilon (1/neuron)');

subplot(2,1,1);
title('[I_0N]/(1+N\epsilon) fit parameters');
%plt_fit2 = Plot();


%%

I_shuffle_mean = fi_m - fi_sh_m;
I_shuffle_errb = sqrt(fi_e.^2 + fi_sh_e.^2);

%I_diag_mean = te_fish_mean - te_fish_mean_preshuf;
%I_diag_errb = sqrt(te_fish_errb.^2 + te_fish_errb_preshuf.^2);

figure;
hold on;

%errorbar(cell_nums, I_diag_mean, I_diag_errb);
errorbar(n_cells, I_shuffle_mean, I_shuffle_errb);

xlabel('Number of cells');
ylabel('\Delta Fisher information (cm^{-2})');
title('Information difference compared to shuffled data');
%legend '\DeltaI_{diag}' '\DeltaI_{shuffle}' 'Location' 'best'
legend '\DeltaI_{shuffle}' 'Location' 'best'

plt_diff = Plot();
%plt_diff.Colors{2} = plt.Colors{3};
%plt_diff.Colors{1} = plt.Colors{2};
plt_diff.XMinorTick = 'off';
plt_diff.YMinorTick = 'off';
plt_diff.ShowBox = 'off';


%%
function [means, errbs, n_cells] = fetch_from_rec(rec, out_type, trainset)
if nargin == 2
    trainset = false;
end

if strcmp(out_type, 'fisher_info')
    func = @(x) 1./x;
    out_type = 'MSE';
else
    func = @(x) x;
end

if trainset
    out_type = [out_type '_shuf'];
end

outp = arrayfun(@(x) func(mean(x.(out_type))), rec);
means = mean(outp, 2);
errbs = std(outp, [], 2)./sqrt(size(rec,2));
n_cells = [rec(:,1).num_cells].';
end