function f2d
%% Figure 2d
% Depends on file: decoder_by_bins2_agg_210909-141200_0.mat
% Raw numbers for:
% - LineX
% - LineRealY, ErrorBarsRealY
% - LineShuffledY, ErrorBarsShuffledY
rng(0);

data = decoders_by_bins;
make_xlsx(data, 'f2d');
end

%% helper functions
function data = decoders_by_bins
load('decoder_by_bins2_agg_210909-141200_0.mat', 'res');

m = squeeze(mean(cat(3, res.me_bins_max), 1)).';
m_s = squeeze(mean(cat(3, res.me_bins_shuf_max), 1)).';
n_sess = size(m, 1);

figure;
errorbar(1:20, mean(m), 1.96*std(m)./sqrt(n_sess), 'b');
hold on;

data.LineX = 1:20;
data.LineRealY = mean(m);
data.ErrorBarsRealY = 1.96*std(m)./sqrt(n_sess);

errorbar(1:20, mean(m_s), 1.96*std(m_s)./sqrt(n_sess), 'r');

data.LineShuffledY = mean(m_s);
data.ErrorBarsShuffledY = 1.96*std(m_s)./sqrt(n_sess);

xlabel 'Spatial bin'
ylabel 'RMS error (cm)'

text(15, 15, 'Real', 'Color', 'b');
text(10, 3.5, 'Shuffled', 'Color', 'r');

figure_format([1.5 1.2]);
end