function f5h
%% Figure 5h
% Depends on files: default_store.mat
% Raw numbers for:
% - Values2021Real, Values2021Shuffled
% - Values2026Real, Values2026Shuffled
rng(0);

data = asnr_boxplot;
make_xlsx(data, 'f5h');
end

%% helper functions
function data = asnr_boxplot
org = Org;
org.init;

mouse_low = 'Mouse2021';
mouse_high = 'Mouse2026';

figure;
subplot(2,2, [2 4]);
get_prop = @(p, m) org.sess_prop.(p)(strcmp(org.mouse, m));
get_mouse = @(m) org.mouse(strcmp(org.mouse, m));

boxplot(log10([get_prop('asymp_snr', mouse_low);...
    get_prop('asymp_snr', mouse_high)]),...
    [get_mouse(mouse_low), get_mouse(mouse_high)],...
    'Colors', 'b', 'Symbol', '+');
hold on;

data.Values2021Real = get_prop('asymp_snr', mouse_low);
data.Values2026Real = get_prop('asymp_snr', mouse_high);

boxplot(log10([get_prop('asymp_snr_shuf', mouse_low);...
    get_prop('asymp_snr_shuf', mouse_high)]),...
    [get_mouse(mouse_low), get_mouse(mouse_high)],...
    'Colors', 'r', 'Symbol', '+');

data.Values2021Shuffled = get_prop('asymp_snr_shuf', mouse_low);
data.Values2026Shuffled = get_prop('asymp_snr_shuf', mouse_high);

ylabel 'log_{10} asymptotic SNR'
text(1.5, 0, 'Real', 'Color', 'b', 'HorizontalAlignment','center');
text(1.5, 1.5, 'Shuffled', 'Color', 'r', 'HorizontalAlignment','center');

drawArrow = @(x,y,varargin) quiver( x(1),y(1),x(2)-x(1),y(2)-y(1),0, varargin{:} );
m_color = DecodeTensor.mcolor({mouse_low});
m_color = m_color{1};
b = mean(log10(get_prop('asymp_snr', mouse_low)));
t = mean(log10(get_prop('asymp_snr_shuf', mouse_low)));
drawArrow([1 1]-0.05,[b t], 'Color', m_color);
text(1 + 0.05, (b+t)/2, sprintf('x%.2g', 10^(t-b)), 'Color', m_color);

m_color = DecodeTensor.mcolor({mouse_high});
m_color = m_color{1};
b = mean(log10(get_prop('asymp_snr', mouse_high)));
t = mean(log10(get_prop('asymp_snr_shuf', mouse_high)));
drawArrow([2 2]-0.05,[b t], 'Color', m_color);
text(2 + 0.1, (b+t)/2, sprintf('x%.2g', 10^(t-b)), 'Color', m_color);

set(gca, 'Box', 'off');
end