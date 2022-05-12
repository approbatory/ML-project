function f4a
%% Figure 4a
% Depends on files: eigen_snr_crossval_save.mat
% Raw numbers for:
% - LineX, LineTrainY, LineTestY
% - ShadedTrainY, ShadedTestY
rng(0);

data = eigen_snr;
make_xlsx(data, 'f4a');
end

%% helper functions
function data = eigen_snr
load eigen_snr_crossval_save.mat dp2_test dp2_train

filt = SessManager.highqual_filt_from_usable;
%mouse_ = org.mouse(filt);
dp2_test = dp2_test(filt);
dp2_train = dp2_train(filt);
%{
noise_test = noise_test(filt);
noise_train = noise_train(filt);
signal_test = signal_test(filt);
signal_train = signal_train(filt);
%}
%%
[dp2_test_agg, dp2_test_agg_sem] = jagged_agg(dp2_test);
[dp2_train_agg, dp2_train_agg_sem] = jagged_agg(dp2_train);

%{
[signal_test_agg, signal_test_agg_sem] = jagged_agg(signal_test);
[signal_train_agg, signal_train_agg_sem] = jagged_agg(signal_train);

[noise_test_agg, noise_test_agg_sem] = jagged_agg(noise_test);
[noise_train_agg, noise_train_agg_sem] = jagged_agg(noise_train);
%}

%mean over random samples
dp2_test_agg_mean = squeeze(mean(dp2_test_agg,1));
dp2_test_agg_mean_sem = squeeze(sqrt(mean(dp2_test_agg_sem.^2,1) + (std(dp2_test_agg,0,1)./sqrt(size(dp2_test_agg,1))).^2));

dp2_train_agg_mean = squeeze(mean(dp2_train_agg,1));
dp2_train_agg_mean_sem = squeeze(sqrt(mean(dp2_train_agg_sem.^2,1) + (std(dp2_train_agg,0,1)./sqrt(size(dp2_train_agg,1))).^2));

%{
signal_test_agg_mean = squeeze(mean(signal_test_agg,1));
signal_test_agg_mean_sem = squeeze(sqrt(mean(signal_test_agg_sem.^2,1) + (std(signal_test_agg,0,1)./sqrt(size(signal_test_agg,1))).^2));

signal_train_agg_mean = squeeze(mean(signal_train_agg,1));
signal_train_agg_mean_sem = squeeze(sqrt(mean(signal_train_agg_sem.^2,1) + (std(signal_train_agg,0,1)./sqrt(size(signal_train_agg,1))).^2));


noise_train_agg_mean = squeeze(mean(noise_train_agg,1));
noise_train_agg_mean_sem = squeeze(sqrt(mean(noise_train_agg_sem.^2,1) + (std(noise_train_agg,0,1)./sqrt(size(noise_train_agg,1))).^2));

noise_test_agg_mean = squeeze(mean(noise_test_agg,1));
noise_test_agg_mean_sem = squeeze(sqrt(mean(noise_test_agg_sem.^2,1) + (std(noise_test_agg,0,1)./sqrt(size(noise_test_agg,1))).^2));
%}
%% dp2
figure;

H_TR = serrorbar(median(dp2_train_agg_mean), median(dp2_train_agg_mean_sem) * 1.96, 'k:');
hold on;

a_ = median(dp2_train_agg_mean);
data.LineX = 1:numel(a_);
data.LineTrainY = a_;
data.ShadedTrainY = median(dp2_train_agg_mean_sem) * 1.96;

H_TE = serrorbar(median(dp2_test_agg_mean), median(dp2_test_agg_mean_sem) * 1.96, 'b');

data.LineTestY = median(dp2_test_agg_mean);
data.ShadedTestY = median(dp2_test_agg_mean_sem) * 1.96;

[~, argmax_pc_snr] = max(median(dp2_test_agg_mean));
med_dp2_test = median(dp2_test_agg_mean);
fprintf('Sum of dp2 from 1 to %d (the max) = %f\n which is %f%% of the total %f\n', argmax_pc_snr, sum(med_dp2_test(1:argmax_pc_snr)), 100 * sum(med_dp2_test(1:argmax_pc_snr))/sum(med_dp2_test), sum(med_dp2_test));
l_ = line([argmax_pc_snr argmax_pc_snr], ylim, 'Color', 'r');

x_ = median(dp2_test_agg_mean);
cum_x_ = cumsum(x_) ./ sum(x_);
ix_most = find(cum_x_ >= 0.95, 1);
l2_ = line([ix_most ix_most], ylim, 'Color', 'g');


legend([l_, l2_],...
    sprintf('Max SNR (PC%d)', argmax_pc_snr),...
    sprintf('95%% SNR (PC1-%d)', ix_most));

legend boxoff
legend location north


xlabel 'PC index'
ylabel 'Eigenmode SNR'
ylim([0 0.11]);
end

function [m, s] = jagged_agg(X) %jagged on dim 3
    N = length(X);
    max_dim = max(cellfun(@(x)size(x,3), X));
    Y = cell(1,max_dim);
    for i = 1:max_dim
        for j = 1:N
            my_max_dim = size(X{j},3);
            if i > my_max_dim
                continue;
            end
            Y{i}{end+1} = X{j}(:,:,i);
        end
    end
    
    for i = 1:max_dim
        C = cat(3, Y{i}{:});
        num_samples = size(C,3);
        m(:,:,i) = mean(C,3);
        s(:,:,i) = std(C,0,3)./sqrt(num_samples);
    end
end