function [dmu_k, dmu_k2, bins] = signal_density_demo(o, save_fig)
if ~exist('save_fig', 'var')
    save_fig = false;
end

dmu_k = cell2mat(o.dmus);
dmu_k2 = dmu_k.^2;

bins = [(1:19), (21:39)];
bin_labels = [cellfun(@(x)[num2str(x) 'R'],(num2cell(1:19)), 'UniformOutput', false),...
    cellfun(@(x)[num2str(x) 'L'],(num2cell(1:19)), 'UniformOutput', false)];

dmu_k = sort_each_col(dmu_k);
dmu_k2 = sort_each_col(dmu_k2);

show_ticks = 1:4:38;

figure('Position', [0.1517    0.4183    1.8607    0.7707]*1e3);
subplot(4,3,[1 4 7 10]);
imagesc(dmu_k);
colormap(bluewhitered);
set(gca, 'XTick', show_ticks);
set(gca, 'XTickLabel', bin_labels(show_ticks));
xlabel 'Spatial bins'
ylabel 'Neurons (sorted)'
h_ = colorbar;
h_.Title.String = '\Delta\mu_k (\DeltaF/F)';

subplot(4,3,[2 5 8 11]);
plot(dmu_k2, '.-');
%set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
xlabel 'Neurons (sorted)'
ylabel '(\Delta\mu_k)^2 (\DeltaF/F)^2'
%imagesc(dmu_k2);
%set(gca, 'XTick', show_ticks);
%set(gca, 'XTickLabel', bin_labels(show_ticks));
%set(gca,'colorscale','log');
%xlabel 'Spatial bins'
%ylabel 'Neurons (sorted)'
%h_ = colorbar;
%h_.Title.String = '(\Delta\mu_k)^2 (\DeltaF/F)^2';

subplot(4,3,[3 6 9]);%change indices
dmu_k2_normed = dmu_k2 ./ sum(dmu_k2);
plot((1:o.num_neurons)./o.num_neurons, dmu_k2_normed, '.',...
    'Color', [1 1 1]*0.6);
hold on
hp = plot((1:o.num_neurons)./o.num_neurons, mean(dmu_k2_normed,2), '-xk');
%set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
xlim([0 1]);
xlabel 'Fraction of neurons (sorted)'
ylabel '(\Delta\mu_k)^2 normalized'
%imagesc(dmu_k2_normed);
%set(gca, 'XTick', show_ticks);
%set(gca, 'XTickLabel', bin_labels(show_ticks));
%set(gca,'colorscale','log');
%xlabel 'Spatial bins'
%ylabel 'Neurons (sorted)'
%h_ = colorbar;
%h_.Title.String = '(\Delta\mu_k)^2 normalized';
ipr = 1 ./ sum(dmu_k2_normed.^2);
sig_dens = ipr ./ o.num_neurons;
hold on;
for i = 1:numel(sig_dens)
    s_= sig_dens(i);
    hl_ = line([s_ s_], ylim, 'Color', [1 1 1]*0.8); 
end
s_ = mean(sig_dens);
hl = line([s_ s_], ylim, 'Color', 'k', 'LineWidth', 2);
legend([hp hl_ hl], 'Mean over spatial bins',...
    'Signal density', 'Mean signal density');
legend box off
legend location northwest


subplot(4,3,12);
plot(bins, sig_dens, '-x');
set(gca, 'XTick', show_ticks);
set(gca, 'XTickLabel', bin_labels(show_ticks));
xlabel 'Spatial bins'
ylabel 'Signal density'
hold on;
l_ = refline([0 mean(ipr ./ size(dmu_k2,1))]);
l_.Color = 'k';
legend 'Density per bin' 'Mean density'
legend boxoff
legend Location best

suptitle([o.dt.mouse_name esc(o.dt.source_path(72:end-33))]);

if save_fig
    fname = sprintf('signal_density_demo_%s.fig', o.dt.mouse_name);
    savefig(fname);
end

return;
figure;
subplot(1,2,1);
hold on;
for j = 1:size(dmu_k,2)
    histogram(dmu_k(:,j));
    [~, p_(j)] = lillietest(dmu_k(:,j));
end
title(sprintf('Showing %d hists', size(dmu_k,2)));
xlabel '\Delta\mu_k (\DeltaF/F)'
ylabel 'Frequency'
subplot(1,2,2);
histogram(p_);
xlabel 'Normality test p-val'
title(sprintf('p = %f +- %f', mean(p_), sem(p_)));

figure;
subplot(1,2,1);
hold on;
for j = 1:size(dmu_k2,2)
    histogram(log(dmu_k2(:,j)));
    [~, p_(j)] = lillietest(log(dmu_k2(:,j)));
end
title(sprintf('Showing %d hists', size(dmu_k,2)));
xlabel 'log \Delta\mu_k^2 log[(\DeltaF/F)^2]'
ylabel 'Frequency'
subplot(1,2,2);
histogram(p_);
xlabel 'Normality test p-val'
title(sprintf('p = %f +- %f', mean(p_), sem(p_)));
end


function X = sort_each_col(X)
[~, cols] = size(X);

for j = 1:cols
    X(:,j) = sort(X(:,j), 'descend');
end
end

function X = sort_by_col1(X)
[~, order] = sort(X(:,1), 'descend');
X = X(order, :);
end