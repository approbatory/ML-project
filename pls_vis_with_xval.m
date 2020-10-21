sm = SessManager;
d = sm.cons_usable(69);
[T1, d1, T2, d2] = DecodeTensor.holdout_half(d.data_tensor, d.tr_dir);
[X1, ks1] = DecodeTensor.tensor2dataset(T1, d1);
[X2, ks2] = DecodeTensor.tensor2dataset(T2, d2);
X2_shuf = shuffle(X2, ks2);

X1z = zscore(X1);
X2z = zscore(X2);
X2z_shuf = zscore(X2_shuf);

[XS1, stats, origin] = Utils.pls_short(X1z, [ceil(ks1/2), mod(ks1,2)]);

apply_proj = @(x) (x - mean(x)) * stats.W;

XS2 = apply_proj(X2z);
origin2 = -mean(X2z) * stats.W;
XS2_shuf = apply_proj(X2z_shuf);
origin2_shuf = -mean(X2z_shuf) * stats.W;

figure;
scatter(XS2(:,1), XS2(:,2), 8, Utils.colorcode(ceil(ks2/2)), 'filled', 'MarkerFaceAlpha', 0.02); hold on;
scatter(origin2(1), origin2(2), 10, 'k');
xlabel PLS1; ylabel PLS2;
axis equal; axis tight;
xl_ = xlim;
yl_ = ylim;
set(gca, 'XTick', []);
set(gca, 'YTick', []);
figure_format('boxsize', [1 1.2], 'fontsize', 6);
print('-dpng', '-r1800', 'figure2_pdf/demo/PLS_adjacent_XVAL.png');

figure;
scatter(XS2_shuf(:,1), XS2_shuf(:,2), 8, Utils.colorcode(ceil(ks2/2)), 'filled', 'MarkerFaceAlpha', 0.02); hold on;
scatter(origin2_shuf(1), origin2_shuf(2), 10, 'k');
xlabel PLS1; ylabel PLS2;
axis equal; axis tight;
xlim(xl_); ylim(yl_);
set(gca, 'XTick', []);
set(gca, 'YTick', []);
figure_format('boxsize', [1 1.2], 'fontsize', 6);
print('-dpng', '-r1800', 'figure2_pdf/demo/PLS_adjacent_shuf_XVAL.png');

