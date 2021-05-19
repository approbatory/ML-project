p = Pub(14, 8, 'rows', 2, 'columns', 3, 'start_char', 'f');

%3f
%{
p.panel(1,...
    'xlab', 'Mouse name (Mouse20XX)',...
    'ylab', 'Asymptotic SNR');
org.boxplot('asymp_snr', 'asymp_snr_shuf', true, false);
%}

%3g
p.panel(2,...
    'xlab', 'Single-cell SNR',...
    'ylab', '{\itI}_0 fit value (cm^{-2}\cdotneurons^{-1})');
org.correlogram('single_dp2', 'I0', true, true);
p.format;
Utils.fix_exponent(gca, 'y', 0);

%3h
p.panel(3,...
    'xlab', 'Asymptotic SNR',...
    'ylab', sprintf('Asymptotic fit IMSE (cm^{-2})'));
org.correlogram('asymp_snr', 'I0N', true, true);

%3i
%{
p.panel(4,...
    'xlab', 'PC index',...
    'ylab', 'Eigenmode SNR');
org.eigen_snr_crossval_aggregation(true);
ylim([0 0.11]);
%}

%3j
%{
p.panel(5,...
    'xlab', 'PC index',...
    'ylab', 'cos^2(PC_i, \Delta\mu)');
org.area_between_cos2(24);
%}

%3k
MAX_DIM = 24;
p.panel(6,...
    'xlab', sprintf('Area between PC_{[1-%d]}', MAX_DIM),...
    'ylab', '1/{\itN} fit value (neurons^{-1})');
name = sprintf('delta_cos2_area_%d', MAX_DIM);
org.correlogram(name, 'invN50', true, true, true);
Nvals = [200 100 65 50];
Nlabels = arrayfun(@(x)['1/' num2str(x)],Nvals,'UniformOutput',false);
set(gca, 'YTick', 1./Nvals);
set(gca, 'YTickLabels', Nlabels);

p.format;
p.print('.', 'f3_update', true);
savefig('f3_update');