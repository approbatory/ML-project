% Pairwise correlated fluctuations between neurons as a function of the distance between the RFs’ peaks (3,941,140 pairs from 110 sessions from 12 mice). There is a positive mean correlation between RFs with peaks separated by up to an ? of the track length.
fig_name = 'f1d_noiseCorrByPeakDistance';

rng default

org = Org;
org.load_definitions;
%%
figure;
ncd_r = cell2mat(org.fetch('nc_by_dist_r')');
r_h = serrorbar(0:19, mean(ncd_r), sem(ncd_r).*1.96, 'lineprops', {'-', 'Color', [0 0.5 0]});
hold on;
ncd_l = cell2mat(org.fetch('nc_by_dist_l')');
l_h = serrorbar(0:19, mean(ncd_l), sem(ncd_l).*1.96, 'lineprops', {'-', 'Color', [0.8 0.5 0]});

l_ = refline(0,0); l_.Color = 'k';
xlabel 'RFs'' peak distance'
ylabel 'Pairwise noise correlations'
legend([r_h.mainLine l_h.mainLine], 'Running right', 'Running left');
legend boxoff

figure_format([1.3 1.3]);
set(gca, 'XTick', [0 10 20]);
set(gca, 'XTickLabel', {'0', 'L/2', 'L'});

render_fig(fig_name);

%% per mouse:

mice = unique(org.mouse);
figure;

for i = 1:numel(mice)
    subplot(3,4,i);
    my_mouse = mice{i};
    ncd_r = permute(org.mouse_all_sess('nc_by_dist_r', my_mouse), [3 2 1]);
    r_h = serrorbar(0:19, mean(ncd_r), sem(ncd_r).*1.96, 'lineprops', {'-', 'Color', [0 0.5 0]});
    hold on;
    
    ncd_l = permute(org.mouse_all_sess('nc_by_dist_l', my_mouse), [3 2 1]);
    l_h = serrorbar(0:19, mean(ncd_l), sem(ncd_l).*1.96, 'lineprops', {'-', 'Color', [0.8 0.5 0]});
    
    l_ = refline(0,0); l_.Color = 'k';
    xlabel 'RFs'' peak distance'
    ylabel 'Pairwise noise correlations'
    legend([r_h.mainLine l_h.mainLine], 'Running right', 'Running left');
    legend boxoff
    
    %figure_format([1.3 1.3]);
    set(gca, 'XTick', [0 10 20]);
    set(gca, 'XTickLabel', {'0', 'L/2', 'L'});
    ylim([-0.1 0.1]);
    
    text(2, -0.06, my_mouse);
    fst = @(x)x(1);
    text(2, -0.07, sprintf('Same bin corr: %.3g (l) , %.3g (r)',...
        fst(mean(ncd_l)), fst(mean(ncd_r))));
    text(2, -0.08, sprintf('Dist to neg. corr: %d/20L (l) , %d/20L (r)',...
        largest_before_dip(mean(ncd_l)), largest_before_dip(mean(ncd_r))));
    
    same_bin_corr_l(i) = fst(mean(ncd_l)); 
    same_bin_corr_l_conf(i) = fst(sem(ncd_l).*1.96);
    same_bin_corr_r(i) = fst(mean(ncd_r)); 
    same_bin_corr_r_conf(i) = fst(sem(ncd_r).*1.96);
    
    bin_dist_to_neg_corr_l(i) = largest_before_dip(mean(ncd_l));
    bin_dist_to_neg_corr_r(i) = largest_before_dip(mean(ncd_r));
end

corr_properties_table = table(mice(:), same_bin_corr_l(:), same_bin_corr_l_conf(:),...
    same_bin_corr_r(:), same_bin_corr_r_conf(:), bin_dist_to_neg_corr_l(:), bin_dist_to_neg_corr_r(:),...
    'VariableNames', {'Mouse', 'SameBinCorrLeft', 'SameBinCorrLeftConf',...
    'SameBinCorrRight', 'SameBinCorrRightConf', 'PerMouseBinDistToNegativeCorrLeft',...
    'PerMouseBinDistToNegativeCorrRight'});

function n = largest_before_dip(x)
    fst_neg = find(x < 0, 1);
    if isempty(fst_neg)
        n = numel(x);
        return;
    end
    n = fst_neg - 1;
end