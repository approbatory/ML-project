d = DecodeTensor.cons_filt(70);

data = d.data_tensor;
data_shuf = DecodeTensor.shuffle_tensor(d.data_tensor, d.tr_dir);

%%
map_both = squeeze(data(:,10,:))';

map_l = squeeze(data(:,10,d.tr_dir==-1))';
map_l_next = squeeze(data(:,18,d.tr_dir==-1))';
map_r = squeeze(data(:,10,d.tr_dir== 1))';
map_r_next = squeeze(data(:,18,d.tr_dir== 1))';

map_both_shuf = squeeze(data_shuf(:,10,:))';
map_l_shuf = squeeze(data_shuf(:,10,d.tr_dir==-1))';
map_r_shuf = squeeze(data_shuf(:,10,d.tr_dir== 1))';

%%
figure;
idx = kmeans(map_both, 2);
subplot(2,1,1);
Utils.pca_plot(map_both, d.tr_dir);
subplot(2,1,2);
Utils.pca_plot(map_both, idx);
missed_trials = sum((d.tr_dir + 1)/2+1 ~= idx'); %only 3 trials were misclassified
%%
data_ = data - mean(data,3);
data_shuf_ = data_shuf - mean(data_shuf,3);

figure;
subplot(2,1,1);
for b_i = 1:20
    Utils.pca_plot( squeeze(data_(:,b_i,:))', d.tr_dir);
    hold on
end
title real

subplot(2,1,2);
for b_i = 1:20
    Utils.pca_plot( squeeze(data_shuf_(:,b_i,:))', d.tr_dir);
    hold on
end
title shuf

colormap jet

%%
figure;
Utils.pca_plot([map_l; map_l_next; map_r; map_r_next], vr([1 2 3 4], [size(map_l,1), size(map_l_next,1), size(map_r,1), size(map_r_next,1)]));
colormap jet

%%
[XX, kk] = DecodeTensor.tensor2dataset(data, d.tr_dir);
Utils.pca_plot(XX, mod(kk,2));
%%
Utils.pca_plot(XX(mod(kk,2)==0,:), kk(mod(kk,2)==0));
%%

[XX, kk] = DecodeTensor.tensor2dataset(data - mean(data,3), d.tr_dir);
Utils.pca_plot(XX, mod(kk,2));
%%
Utils.pca_plot(XX(mod(kk,2)==0,:), kk(mod(kk,2)==0));
%%
d = DecodeTensor.cons_filt(70);
parfor rep_i = 1:20
    [sub_data, d1] = DecodeTensor.holdout_half(d.data_tensor, d.tr_dir);
    sub_data_shuf = DecodeTensor.shuffle_tensor(sub_data, d1);
    %sub_data = d.data_tensor;
    [XX_, kk] = DecodeTensor.tensor2dataset(sub_data, d1);
    [XX_s, kk_s] = DecodeTensor.tensor2dataset(sub_data_shuf, d1);
    for n_cl = 1:10
        [~, ~, sumd{rep_i, n_cl}] = kmeans(XX_, n_cl, 'Distance', 'cosine');
        num_samples(rep_i) = size(XX_,1);
        [~, ~, sumd_s{rep_i, n_cl}] = kmeans(XX_s, n_cl, 'Distance', 'cosine');
    end
end
%
cos_err = 1-cellfun(@sum, sumd)./num_samples';
cos_err_shuf = 1-cellfun(@sum, sumd_s)./num_samples';
%cos_err = cellfun(@sum, sumd)./num_samples';
%cos_err_shuf = cellfun(@sum, sumd_s)./num_samples';
%%

figure;
shadedErrorBar(1:10, mean(cos_err_shuf), sem(cos_err_shuf), 'lineprops', 'r'); hold on
shadedErrorBar(1:10, mean(cos_err), sem(cos_err), 'lineprops', 'b'); hold on;

xlabel 'Number of clusters, {\itk}'
ylabel(sprintf('Mean cosine to\ncluster centroid'))
line([2 2], [0 1], 'LineStyle', ':', 'Color', 'k');
ylim([0.25 0.7]);
xlim([1 10]);

text(6, 0.45, 'Real', 'Color', 'b');
text(6, 0.35, 'Shuffled', 'Color', 'r');
figure_format;
Utils.printto('supplements_pdf', 'kmeans_test.pdf');
%TODO: run on all sessions!
%%

function res = vr(v, r)
res = [];
assert(numel(v)==numel(r));
for i = 1:numel(v)
    res = [res, repmat(v(i), 1, r(i))];
end
end