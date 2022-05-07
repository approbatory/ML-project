function f1b
%% Figure 1b
% Raw numbers for:
% - LineX, LineY, ErrorBarsY
% - ScatterPointsX, ScatterPointsY
rng(0);

data = plot_example_pf('Mouse2022', 326, 4);
make_xlsx(data, 'f1b');
end

%% helper functions
function data = plot_example_pf(mouse, ix, use_jitter)
if ~exist('use_jitter', 'var')
    use_jitter = false;
end
d = SessManager.load_special(mouse, 'events_transients');
S = load(d.source_path);

%ix = 326;
ntype = 'events_transients';
tr = S.tracesEvents.(ntype)(:,ix);

pos = S.tracesEvents.position(:,1);
p_low = prctile(pos,1);
p_high = prctile(pos,99);
pos(pos < p_low) = p_low;
pos(pos > p_high) = p_high;
pos = (pos - p_low) ./ (p_high - p_low) .* 120;
v_filt = diff(pos) > 4/20;

tr = tr(v_filt);
pos = pos(v_filt);
[i,j,v] = find(tr);

p = pos(i);
jitter = use_jitter.*120*0.01*(rand(size(p))-0.5);
figure;
scatter(p+jitter,v, 5, 'b', 'filled');  %%This line plots the blue points

data.ScatterPointsX = p+jitter;
data.ScatterPointsY = v;

alpha(0.1);
xlim([0 120]);
set(gca, 'XTick', [0 60 120]);
set(gca, 'YTick', [0,2,4,6,8]);
set(gca, 'Box', 'off');
xlabel 'Track position (cm)'
ylabel 'Scaled events'
title(sprintf('Example neuron (%d rightward trials)', sum(d.tr_dir == 1)));

mean_act = mean(d.data_tensor(ix,:,d.tr_dir==1),3);
std_act = std(d.data_tensor(ix,:,d.tr_dir==1),0,3);
bin_centers = 3:6:120;
hold on;
errorbar(bin_centers, mean_act, std_act, 'b');

data.LineX = bin_centers;
data.LineY = mean_act;
data.ErrorBarsY = std_act;

fprintf('The standard deviation of activity per bin is: %f +- %f (std)\n', mean(std_act), std(std_act));

figure_format([1.8 1.8]);
end