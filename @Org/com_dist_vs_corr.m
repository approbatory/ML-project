function com_dist_vs_corr(o)

agg = @(x) cell2mat(Utils.cf_(@get_pairs, o.fetch(x)'));

total_noise_corr_right = agg('noise_corr_avg_right');
total_noise_corr_left = agg('noise_corr_avg_left');

total_sig_corr_right = agg('sig_corr_right');
total_sig_corr_left = agg('sig_corr_left');

total_corr_right = agg('tot_corr_right');
total_corr_left = agg('tot_corr_left');

total_com_dist_right = agg('com_dist_right');
total_com_dist_left = agg('com_dist_left');

total_mode_dist_right = agg('mode_dist_right');
total_mode_dist_left = agg('mode_dist_left');

bin_params = {0,10,40};

figure;
subplot(3,2,1);
[d,c,cs] = show_vars(total_com_dist_right, total_noise_corr_right,...
    bin_params{:});
serrorbar(d,c,cs);
l_ = refline(0,0); l_.Color = 'k';
xlabel 'COM dist. right (bins)'
ylabel 'Noise corr. right'

subplot(3,2,2);
[d,c,cs] = show_vars(total_com_dist_left, total_noise_corr_left,...
    bin_params{:});
serrorbar(d,c,cs);
l_ = refline(0,0); l_.Color = 'k';
xlabel 'COM dist. left (bins)'
ylabel 'Noise corr. left'

subplot(3,2,3);
[d,c,cs] = show_vars(total_com_dist_right, total_sig_corr_right,...
    bin_params{:});
serrorbar(d,c,cs);
l_ = refline(0,0); l_.Color = 'k';
xlabel 'COM dist. right (bins)'
ylabel 'Signal corr. right'

subplot(3,2,4);
[d,c,cs] = show_vars(total_com_dist_left, total_sig_corr_left,...
    bin_params{:});
serrorbar(d,c,cs);
l_ = refline(0,0); l_.Color = 'k';
xlabel 'COM dist. left (bins)'
ylabel 'Signal corr. left'

subplot(3,2,5);
[d,c,cs] = show_vars(total_com_dist_right, total_corr_right,...
    bin_params{:});
serrorbar(d,c,cs);
l_ = refline(0,0); l_.Color = 'k';
xlabel 'COM dist. right (bins)'
ylabel 'Total corr. right'

subplot(3,2,6);
[d,c,cs] = show_vars(total_com_dist_left, total_corr_left,...
    bin_params{:});
serrorbar(d,c,cs);
l_ = refline(0,0); l_.Color = 'k';
xlabel 'COM dist. left (bins)'
ylabel 'Total corr. left'

suptitle('All sessions, using COM distance');



bin_params = {0,20,20};
figure;
subplot(3,2,1);
[d,c,cs] = show_vars(total_mode_dist_right, total_noise_corr_right,...
    bin_params{:});
serrorbar(d,c,cs);
l_ = refline(0,0); l_.Color = 'k';
xlabel 'Peak dist. right (bins)'
ylabel 'Noise corr. right'

subplot(3,2,2);
[d,c,cs] = show_vars(total_mode_dist_left, total_noise_corr_left,...
    bin_params{:});
serrorbar(d,c,cs);
l_ = refline(0,0); l_.Color = 'k';
xlabel 'Peak dist. left (bins)'
ylabel 'Noise corr. left'

subplot(3,2,3);
[d,c,cs] = show_vars(total_mode_dist_right, total_sig_corr_right,...
    bin_params{:});
serrorbar(d,c,cs);
l_ = refline(0,0); l_.Color = 'k';
xlabel 'Peak dist. right (bins)'
ylabel 'Signal corr. right'

subplot(3,2,4);
[d,c,cs] = show_vars(total_mode_dist_left, total_sig_corr_left,...
    bin_params{:});
serrorbar(d,c,cs);
l_ = refline(0,0); l_.Color = 'k';
xlabel 'Peak dist. left (bins)'
ylabel 'Signal corr. left'

subplot(3,2,5);
[d,c,cs] = show_vars(total_mode_dist_right, total_corr_right,...
    bin_params{:});
serrorbar(d,c,cs);
l_ = refline(0,0); l_.Color = 'k';
xlabel 'Peak dist. right (bins)'
ylabel 'Total corr. right'

subplot(3,2,6);
[d,c,cs] = show_vars(total_mode_dist_left, total_corr_left,...
    bin_params{:});
serrorbar(d,c,cs);
l_ = refline(0,0); l_.Color = 'k';
xlabel 'Peak dist. left (bins)'
ylabel 'Total corr. left'


suptitle('All sessions, using peak (mode) distance');
end


function [edge_mids, mean_vals, sems] = show_vars(X, Y, low, high, count)
    %[~, edges] = histcounts(X);
    edges = linspace(low, high, count + 1);
    [mean_vals, sems] = avg_bins(X,Y,edges);
    edge_mids = (edges(1:end-1) + edges(2:end))/2;
end


function [mean_vals, sems] = avg_bins(X, Y, Xbins)
[X, ord] = sort(X);
Y = Y(ord);
max_X = max(X);

n_bins = numel(Xbins)-1;
[mean_vals, sems] = deal(zeros(1,n_bins));

for i = 1:n_bins
    bin_start = Xbins(i);
    bin_end = Xbins(i+1);
    
    ind_start = find(X >= bin_start, 1);
    assert(~isempty(ind_start), 'could not find start');
    if bin_end >= max_X
        ind_end = numel(X);
    else
        ind_end = find(X > bin_end, 1) - 1;
        assert(~isempty(ind_end), 'could not find end');
    end
    
    selection = Y(ind_start:ind_end);
    mean_vals(i) = mean(selection);
    sems(i) = sem(selection);
end
end