function f5g
%% Figure 5g
% Depends on files: default_store.mat
% Raw numbers for:
% - LineLeftX, LineLeft2021Y, ShadedLeft2021Y, LineLeft2026Y, ShadedLeft2026Y
% - LineRightX, LineRight2021Y, ShadedRight2021Y, LineRight2026Y, ShadedRight2026Y
% - Bar2021X, Bar2026X, Bar2021Y, Bar2026Y
rng(0);

data = nsv_by_bin;
make_xlsx(data, 'f5g');
end

%% helper functions
function data = nsv_by_bin
org = Org;
org.init;

var_name = 'signal_var_norm';
var_label = 'Normalized signal variance';
mouse_low = 'Mouse2021';
mouse_high = 'Mouse2026';
yl_ = [0 6e-4];

figure;

subplot(2,2,1);
mouse_name = mouse_low;
[m,s]=org.mouse_by_bins(var_name, {mouse_name, 'restrict'});
mean_nsv = nanmean(m);
serrorbar(1:20, m(1:20),s(1:20)*1.96);

data.LineRightX = 1:20;
data.LineRight2021Y = m(1:20);
data.ShadedRight2021Y = s(1:20)*1.96;

ax=gca;
ax.XTick = [10 20];
ax.XTickLabel={'60','120'};
ylim(yl_);
Utils.fix_exponent(gca, 'y', 0);
ylabel(var_label);
m_color = DecodeTensor.mcolor({mouse_name});
m_color = m_color{1};
text(5, 4e-4, mouse_name, 'Color', m_color);
l_ = refline(0, mean_nsv);
l_.Color = [0.3 0.3 0.3];
l_.LineStyle = '--';
xlabel 'Running right (cm)'

subplot(2,2,2);
serrorbar(1:20, m(21:40),s(21:40)*1.96);

data.LineLeftX = 1:20;
data.LineLeft2021Y = m(21:40);
data.ShadedLeft2021Y = s(21:40)*1.96;

l_ = refline(0, mean_nsv);
l_.Color = [0.3 0.3 0.3];
l_.LineStyle = '--';
ylim(yl_);
ax=gca;
ax.XTick=[10 20];
ax.XTickLabel={'60','120'};
ax.YAxis.Visible = 'off';
xlabel 'Running left (cm)'



subplot(2,2,3);
mouse_name = mouse_high;
[m,s]=org.mouse_by_bins(var_name, {mouse_name, 'restrict'});
mean_nsv = nanmean(m);
serrorbar(1:20, m(1:20),s(1:20)*1.96);

data.LineRight2026Y = m(1:20);
data.ShadedRight2026Y = s(1:20)*1.96;

ax = gca;
ax.XTick = [10 20];
ax.XTickLabel = {'60', '120'};
ylim(yl_);
Utils.fix_exponent(gca, 'y', 0);
ylabel(var_label);
m_color = DecodeTensor.mcolor({mouse_name});
m_color = m_color{1};
text(5, 1e-4, mouse_name, 'Color', m_color);
l_ = refline(0, mean_nsv);
l_.Color = [0.3 0.3 0.3];
l_.LineStyle = '--';
xlabel 'Running right (cm)'

subplot(2,2,4);
serrorbar(1:20, m(21:40),s(21:40)*1.96);

data.LineLeft2026Y = m(21:40);
data.ShadedLeft2026Y = s(21:40)*1.96;

l_ = refline(0, mean_nsv);
l_.Color = [0.3 0.3 0.3];
l_.LineStyle = '--';
ylim(yl_);
ax=gca;
ax.XTick=[10 20];
ax.XTickLabel={'60','120'};
ax.YAxis.Visible = 'off';
xlabel 'Running left (cm)'

%%
figure;
x_ = org.mouse_all_sess(var_name, {mouse_low, 'restrict'});
x_ = x_(:);
x_(isnan(x_)) = [];

y_ = org.mouse_all_sess(var_name, {mouse_high, 'restrict'});
y_ = y_(:);
y_(isnan(y_)) = [];

bin_width = 5e-5;
subplot(2,2,1);

m_color = DecodeTensor.mcolor({mouse_low});
m_color = m_color{1};
text(0.1*500, 0.2, mouse_low, 'Color', m_color); hold on;
histogram(x_, 'BinWidth', bin_width, 'FaceColor', m_color, 'Normalization', 'probability');

[data.Bar2021Y, data.Bar2021X] = histcounts(x_, 'BinWidth', bin_width, 'Normalization', 'probability');

xlabel(var_label);
ylabel 'Probability'
axis square;
ax1 = gca;
ax1.Box = 'off';
%ax1.XTick = [];
xlim([0 1e-3]);

subplot(2,2,3);
m_color = DecodeTensor.mcolor({mouse_high});
m_color = m_color{1};
text(0.1*500, 0.2, mouse_high, 'Color', m_color); hold on;
histogram(y_, 'BinWidth', bin_width, 'FaceColor', m_color, 'Normalization', 'probability');

[data.Bar2026Y, data.Bar2026X] = histcounts(y_, 'BinWidth', bin_width, 'Normalization', 'probability');

xlabel(var_label);
ylabel 'Probability'
axis square;
ax2 = gca;
ax2.Box = 'off';
linkaxes([ax1 ax2]);
xlim([0 1e-3]);

Utils.fix_exponent(subplot(2,2,1), 'x', 0);
Utils.fix_exponent(subplot(2,2,3), 'x', 0);
end