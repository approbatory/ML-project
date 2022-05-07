function f2g
%% Figure 2g
% Depends on file: decoder_by_bins2_agg_210909-141200_0.mat
% Raw numbers for:
% - LineX
% - LineRealY, LineShuffledY
% - ShadedRealY, ShadedShuffledY
rng(0);

data = track_ends_decoding;
make_xlsx(data, 'f2g');
end

%% helper functions
function data = track_ends_decoding
load decoder_by_bins2_agg_210909-141200_0.mat res

my_mouse = 'Mouse2022';
s_idx = SessManager.special_sessions_usable_index(my_mouse);

bin_groups = {[1 20], [2 3 18 19], [4 5 16 17], [6 7 14 15], [8 9 12 13], [10 11]};
bg_names = {'far end', 'end', 'far side', 'mid side', 'near side', 'center'};

bg_i = 2;

bg = bin_groups{bg_i};

imse_bins = 1 ./ res(s_idx).me_bins.^2;
imse_bins_shuf = 1./ res(s_idx).me_bins_shuf.^2;
n_size = res(s_idx).n_size;

x = squeeze(mean(imse_bins(:, :, bg),3)).';
x_shuf = squeeze(mean(imse_bins_shuf(:, :, bg),3)).';

figure;
shadedErrorBar(n_size, x, {@mean, @(x)sem(x)*1.96}, 'lineprops', 'b');
hold on;

data.LineX = n_size;
data.LineRealY = mean(x);
data.ShadedRealY = sem(x)*1.96;

shadedErrorBar(n_size, x_shuf, {@mean, @(x)sem(x)*1.96}, 'lineprops', 'r');

data.LineShuffledY = mean(x_shuf);
data.LineShuffledY(data.LineShuffledY > 100) = nan;
data.ShadedShuffledY = sem(x_shuf)*1.96;

ylim([0 Inf]);
xlim([0 500]);
xlabel 'Num. cells'
ylabel '1/MSE (cm^{-2})'
title(sprintf('%s from %s', bg_names{bg_i}, my_mouse));
end