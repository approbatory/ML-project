load nbins_agg_200930-174641_0.mat

%%

nbins = res(1).nbins;
N_usable = numel(res);

for i = 1:N_usable
   mean_cerr(i,:) = squeeze(mean(res(i).c_err , 3));
end

%%
figure;

grand_mean_cerr = mean(mean_cerr);
mean_cerr_conf = sem(mean_cerr).*1.96;
serrorbar(nbins, grand_mean_cerr, mean_cerr_conf);
xlabel '# of bins'
ylabel 'RMS decoding error'
title(sprintf('Decoding error by # of bins, 80 samples per session,\nn = %d sessions (95%% CI over sessions)', N_usable));
%%
[lowest_err, i] = min(grand_mean_cerr);
lowest_err_bins = nbins(i);

fprintf('Minimum at %d bins:\n%.2f cm\n', lowest_err_bins, lowest_err);

p = signrank(mean_cerr(:,i-1), mean_cerr(:,i));
p_rs = ranksum(mean_cerr(:,i-1), mean_cerr(:,i));
fprintf('At %d bins, p = %.3e signrank, p = %.3f ranksum test over sessions compared to %d bins\n%.2f cm error\n', nbins(i-1), p, p_rs, nbins(i), grand_mean_cerr(i-1));