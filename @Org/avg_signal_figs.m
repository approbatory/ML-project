%figure;
p = Pub(14, 8, 'rows', 2, 'columns', 3);

%subplot(2,3,1); 
p.panel(1, 'xlab', 'Signal density', 'ylab', 'Asymp. ratio', 'titl', 'Per mouse');
org.correlogram('signal_density', 'asymp_ratio', true, true); title 'Per mouse'; axis square

%subplot(2,3,2); 
p.panel(2, 'xlab', 'Avg. signal', 'ylab', 'Asymp. ratio', 'titl', 'Per mouse');
org.correlogram('avg_signal', 'asymp_ratio', true, true); title 'Per mouse'; axis square

%subplot(2,3,3); 
p.panel(3, 'xlab', 'Avg. signal', 'ylab', 'Signal density', 'titl', 'Per mouse');
org.correlogram('avg_signal', 'signal_density', true, true); title 'Per mouse'; axis square

%subplot(2,3,4); 
p.panel(4, 'xlab', 'Signal density', 'ylab', 'Asymp. ratio', 'titl', 'Per session');
org.correlogram('signal_density', 'asymp_ratio', false, true); title 'Per session'; axis square

%subplot(2,3,5); 
p.panel(5, 'xlab', 'Avg. signal', 'ylab', 'Asymp. ratio', 'titl', 'Per session');
org.correlogram('avg_signal', 'asymp_ratio', false, true); title 'Per session'; axis square

%subplot(2,3,6); 
p.panel(6, 'xlab', 'Avg. signal', 'ylab', 'Signal density', 'titl', 'Per session');
org.correlogram('avg_signal', 'signal_density', false, true); title 'Per session'; axis square

p.format;

%%

p.print('.', 'avg_signal_asymp_ratio');


%%

%figure;
p = Pub(14, 8, 'rows', 2, 'columns', 3);

%subplot(2,3,1); 
p.panel(1, 'xlab', 'Signal density', 'ylab', '{\itI}_0{\itN} ratio', 'titl', 'Per mouse');
org.correlogram('signal_density', 'I0N_ratio', true, true); title 'Per mouse'; axis square

%subplot(2,3,2); 
p.panel(2, 'xlab', 'Avg. signal', 'ylab', '{\itI}_0{\itN} ratio', 'titl', 'Per mouse');
org.correlogram('avg_signal', 'I0N_ratio', true, true); title 'Per mouse'; axis square

%subplot(2,3,3); 
p.panel(3, 'xlab', 'Avg. signal', 'ylab', 'Signal density', 'titl', 'Per mouse');
org.correlogram('avg_signal', 'signal_density', true, true); title 'Per mouse'; axis square

%subplot(2,3,4); 
p.panel(4, 'xlab', 'Signal density', 'ylab', '{\itI}_0{\itN} ratio', 'titl', 'Per session');
org.correlogram('signal_density', 'I0N_ratio', false, true); title 'Per session'; axis square

%subplot(2,3,5); 
p.panel(5, 'xlab', 'Avg. signal', 'ylab', '{\itI}_0{\itN} ratio', 'titl', 'Per session');
org.correlogram('avg_signal', 'I0N_ratio', false, true); title 'Per session'; axis square

%subplot(2,3,6); 
p.panel(6, 'xlab', 'Avg. signal', 'ylab', 'Signal density', 'titl', 'Per session');
org.correlogram('avg_signal', 'signal_density', false, true); title 'Per session'; axis square

p.format;

%%

p.print('.', 'avg_signal_I0N_ratio');