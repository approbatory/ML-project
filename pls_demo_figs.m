load pls_check_output.mat

dp2_train = expand(dp2_train);
dp2_test = expand(dp2_test);

figure;
subplot(2,1,1);
imagesc(squeeze(mean(dp2_train)));
hcb = colorbar;
hcb.Title.String = '(d'')^2';
set(gca, 'ColorScale', 'log');
title 'Train'
xlabel 'PLS dim'
ylabel 'Spatial bin'

subplot(2,1,2);
imagesc(squeeze(mean(dp2_test)));
hcb = colorbar;
hcb.Title.String = '(d'')^2';
set(gca, 'ColorScale', 'log');
title 'Test'
xlabel 'PLS dim'
ylabel 'Spatial bin'

tr_ = squeeze(nanmedian(dp2_train,2));
te_ = squeeze(nanmedian(dp2_test,2));
figure;
errorbar(mean(tr_), sem(tr_), 'k:');
hold on;
errorbar(mean(te_), sem(te_), 'k-');

title 'Median (d'')^2 over bins, unique PLS for adjacent bins'
xlabel 'PLS dim'
ylabel '(d'')^2'

mte_ = mean(te_);
idx = find(diff(mte_) < 0, 1);
line([idx idx], ylim, 'Color', 'k');
legend('Train', 'Test', sprintf('Best dim (%d)', idx));
legend box off
legend location best


function y = expand(x)
y = nan(size(x,1), size(x,2)+2, size(x,3));

y(:,1:19,:) = x(:,1:19,:);
y(:,21:39,:) = x(:,20:38,:);
end