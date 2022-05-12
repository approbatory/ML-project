function f4b
%% Figure 4b
% Depends on files: default_store.mat
% Raw numbers for:
% - LineX, LineRealY, LineShuffledY
% - ShadedRealY, ShadedShuffledY
rng(0);

data = area_between;
make_xlsx(data, 'f4b');
end

%% helper functions
function data = area_between
org = Org;
org.load_definitions;

MAX_DIM = 30;

[l2, l2_sem] = org.all_med_bins('loadings2', 'restrict');
[l2_s, l2_s_sem] = org.all_med_bins('loadings_shuf2', 'restrict');

serrorbar(l2, l2_sem*1.96, 'b');
hold on;

data.LineX = 1:numel(l2);
data.LineRealY = l2;
data.ShadedRealY = l2_sem*1.96;

serrorbar(l2_s, l2_s_sem*1.96, 'r');

data.LineShuffledY = l2_s;
data.ShadedShuffledY = l2_s_sem*1.96;

ylim([0 Inf]);
xlim([1 50]);
patch([1:MAX_DIM, MAX_DIM:-1:1],...
    [l2(1:MAX_DIM)', l2_s(MAX_DIM:-1:1)'],...
    [1 1 1]*0.8, 'FaceAlpha', 0.6, 'EdgeColor', 'none');

xlabel 'PC index'
ylabel 'cos^2(PC_i, \Delta\mu)'


text(20,0.045/3,'Real data', 'Color', 'b');
text(20,0.035/3,'Shuffled', 'Color', 'r');
text(20, 0.025/3,...
    sprintf('Area between PC_{[1-%d]}', MAX_DIM), 'Color', [1 1 1]*0.4);

end