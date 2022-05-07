function f2a
%% Figure 2a
% Raw numbers for:
% - LineX, LineY
% - isExcluded
rng(0);

data = plot_trajectory;
make_xlsx(data, 'f2a');
end

%% helper functions
function data = plot_trajectory
data_source = SessManager.load_special('Mouse2022');
load(data_source.source_path);
pos_raw = tracesEvents.position(:,1);
time = (1:numel(pos_raw))/20;
pos_norm = (pos_raw - prctile(pos_raw, 5))./(prctile(pos_raw, 95) - prctile(pos_raw, 5))*120; %cm
pos_norm(pos_norm < 0) = 0;
pos_norm(pos_norm > 120) = 120;

figure;

plot(time, pos_norm, 'k');

data.LineX = time;
data.LineY = pos_norm;

opt = DecodeTensor.default_opt; opt.total_length = 120; opt.bin_width = 6;
[~, ~, tr_s, tr_e, ~, ~, ~] = DecodeTensor.new_sel(pos_raw, opt);
bad_y = pos_norm;
for tr_i = 1:numel(tr_s)
    bad_y(tr_s(tr_i):tr_e(tr_i)) = nan;
end
hold on;
h = plot(time, bad_y, 'r');

data.isExcluded = ~isnan(bad_y);

legend(h, 'Excluded Frames'); legend boxoff; legend location northoutside

xlim([0 79]);

data_filter = (data.LineX >= 0) & (data.LineX <= 79);
data.LineX = data.LineX(data_filter);
data.LineY = data.LineY(data_filter);
data.isExcluded = data.isExcluded(data_filter);

ylabel 'Position (cm)'
xlabel 'Time (s)'

ylim([0 120]);

figure_format([1.5 1.5]);
end