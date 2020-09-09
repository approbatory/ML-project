function sigdens_plots(org)

figure;

subplot(2,1,1);
mouse_name = 'Mouse2022';
[m,s]=org.mouse_by_bins('signal_density', {mouse_name, 'restrict'});
mean_sig_dens = nanmean(m);
serrorbar(1:20, m(1:20),s(1:20)*1.96);
hold on;
serrorbar(21:40, m(21:40),s(21:40)*1.96);
l_ = refline(0, mean_sig_dens);
l_.Color = [0.3 0.3 0.3];
l_.LineStyle = '--';
ylim([0 0.3]);
ax=gca;
ax.XTick=10:10:40;
ax.XTickLabel={'L/2','L','L/2','L'};
m_color = DecodeTensor.mcolor({mouse_name});
m_color = m_color{1};
text(5, 0.27, mouse_name, 'Color', m_color);

subplot(2,1,2);
mouse_name = 'Mouse2024';
[m,s]=org.mouse_by_bins('signal_density', {mouse_name, 'restrict'});
mean_sig_dens = nanmean(m);
serrorbar(1:20, m(1:20),s(1:20)*1.96);
hold on;
serrorbar(21:40, m(21:40),s(21:40)*1.96);
l_ = refline(0, mean_sig_dens);
l_.Color = [0.3 0.3 0.3];
l_.LineStyle = '--';
ylim([0 0.3]);
ax=gca;
ax.XTick=10:10:40;
ax.XTickLabel={'L/2','L','L/2','L'};
m_color = DecodeTensor.mcolor({mouse_name});
m_color = m_color{1};
text(5, 0.2, mouse_name, 'Color', m_color);

%%
figure;
x_ = org.mouse_all_sess('signal_density', {'Mouse2022', 'restrict'});
x_ = x_(:);
x_(isnan(x_)) = [];

y_ = org.mouse_all_sess('signal_density', {'Mouse2024', 'restrict'});
y_ = y_(:);
y_(isnan(y_)) = [];

bin_width = 0.01;
histogram(x_, 'BinWidth', bin_width, 'FaceColor', 'r', 'Normalization', 'probability');
hold on;
histogram(y_, 'BinWidth', bin_width, 'FaceColor', 'y', 'Normalization', 'probability');
xlabel 'Signal density'
ylabel 'Probability'
%legend 'Mouse2022' 'Mouse2024'

% m_color = DecodeTensor.mcolor({'Mouse2022'});
% m_color = m_color{1};
% text(0.2, 0.23, 'Mouse2022', 'Color', m_color);
% m_color = DecodeTensor.mcolor({'Mouse2024'});
% m_color = m_color{1};
% text(0.2, 0.18, 'Mouse2024', 'Color', m_color);

axis square
figure_format([4 3]/2.5);