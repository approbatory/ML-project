function area_between_cos2(org, MAX_DIM, lean)
if ~exist('lean', 'var')
    lean = true;
end
if ~exist('MAX_DIM', 'var')
    MAX_DIM = 10;
end

[l2, l2_sem] = org.all_med_bins('loadings2', 'restrict');
[l2_s, l2_s_sem] = org.all_med_bins('loadings_shuf2', 'restrict');

serrorbar(l2, l2_sem*1.96, 'b');
hold on;
serrorbar(l2_s, l2_s_sem*1.96, 'r');
ylim([0 Inf]);
xlim([1 50]);
patch([1:MAX_DIM, MAX_DIM:-1:1],...
    [l2(1:MAX_DIM)', l2_s(MAX_DIM:-1:1)'],...
    [1 1 1]*0.8, 'FaceAlpha', 0.6, 'EdgeColor', 'none');

xlabel 'PC index'
ylabel 'cos^2(PC_i, \Delta\mu)'


text(20,0.045,'Real data', 'Color', 'b');
text(20,0.035,'Shuffled', 'Color', 'r');
text(20, 0.025,...
    sprintf('Area between PC_{[1-%d]}', MAX_DIM), 'Color', [1 1 1]*0.4);

%figure_format([4 3]/2.5);



%%
if ~lean
    figure;
    m_list = unique(org.sess_prop.Mouse(org.default_filt(true)));
    for i = 1:numel(m_list)
        subplot(3,4,i);
        [l2, l2_sem] = org.mouse_med_bins('loadings2', {m_list{i},'restrict'});
        [l2_s, l2_s_sem] = org.mouse_med_bins('loadings_shuf2', {m_list{i},'restrict'});
        serrorbar(l2, l2_sem*1.96, 'b');
        hold on;
        serrorbar(l2_s, l2_s_sem*1.96, 'r');
        ylim([0 Inf]);
        xlim([1 50]);
        patch([1:MAX_DIM, MAX_DIM:-1:1],...
            [l2(1:MAX_DIM)', l2_s(MAX_DIM:-1:1)'],...
            [1 1 1]*0.8, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
        
        xlabel 'PC index'
        ylabel 'cos^2(PC_i, \Delta\mu)'
        title(m_list{i});
    end
    
    
    figure;
    org.correlogram(...
        sprintf('delta_cos2_area_%d', MAX_DIM), 'invN50', true, true);
    xlabel(sprintf('\\Delta signal-noise cos^2 overlap (PC1-%d)', MAX_DIM));
    ylabel '1/{\itN} fit value'
    Nvals = [200 100 65 50];
    Nlabels = arrayfun(@(x)['1/' num2str(x)],Nvals,'UniformOutput',false);
    set(gca, 'YTick', 1./Nvals);
    set(gca, 'YTickLabels', Nlabels);
    figure_format([4 3]/2.5);
end