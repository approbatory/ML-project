function [signal_ipr, signal_ipr_sem,...
    spectrum_ipr, spectrum_ipr_sem,...
    asymp_ratio_mean, asymp_ratio_sem, mouse_names] = ...
    all_mice_ipr


[~,mouse_names] = DecodeTensor.special_sess_id_list;
mouse_names = unique(mouse_names);


[signal_ipr, signal_ipr_sem,...
    spectrum_ipr, spectrum_ipr_sem,...
    asymp_ratio_mean, asymp_ratio_sem] = ...
    deal(zeros(1, numel(mouse_names)));


progressbar('mice'); 
for i = 1:numel(mouse_names)
    [sig_ipr, sig_ipr_sem,...
        spec_ipr, spec_ipr_sem,...
        asrat, asrat_sem] = ...
        Cloud.ipr_analysis_for_mouse(mouse_names{i});
    progressbar(i/12);
    
    n = sum(~isnan(sig_ipr));
    signal_ipr(i) = nansum(sig_ipr) ./ n;
    signal_ipr_sem(i) = sqrt(nansum(sig_ipr_sem.^2)) ./ n;
    
    n = numel(spec_ipr);
    spectrum_ipr(i) = sum(spec_ipr) ./ n;
    spectrum_ipr_sem(i) = sqrt(sum(spec_ipr_sem.^2)) ./n;
    
    asymp_ratio_mean(i) = asrat;
    asymp_ratio_sem(i) = asrat_sem;
end

figure;
subplot(1,2,1);
PanelGenerator.plot_regress(signal_ipr, asymp_ratio_mean, ...
    signal_ipr_sem, asymp_ratio_sem, mouse_names, 'k');
title 'Signal IPR vs. asymp ratio'
subplot(1,2,2);
PanelGenerator.plot_regress(spectrum_ipr, asymp_ratio_mean, ...
    spectrum_ipr_sem, asymp_ratio_sem, mouse_names, 'k');
title 'Variance IPR vs. asymp ratio'