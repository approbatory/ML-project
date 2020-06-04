function ipr_analysis_sess

[~, m_] = DecodeTensor.filt_sess_id_list;


load asymp_snr_values.mat
assert(isequal(m_, mouse_names));
[asymp_ratio, asymp_ratio_conf] = ...
    uncertain_divide(asymp_snr, asymp_snr_conf,...
    asymp_snr_shuf, asymp_snr_shuf_conf);

%IPR_savefile = 'corr_KL_ipr_savefile.mat';
IPR_savefile = 'KL_ipr_savefile.mat';
if ~exist(IPR_savefile, 'file')
    
    progressbar('sessions...');
    for idx = 1:numel(m_)
        C = Cloud(idx);
        [spectrum{idx},...
            signal_variances{idx},...
            med_loadings{idx},...
            med_loadings_s{idx},...
            signal_ipr{idx},...
            signal_ipr_s{idx},...
            spectrum_ipr{idx},...
            spectrum_ipr_s{idx},...
            signal_KL{idx},...
            signal_perplexity{idx},...
            area_cos{idx},...
            area_cos2{idx},...
            dmu_ipr{idx},...
            num_neurons(idx),...
            sig_dens{idx},...
            sig_dens_unbiased{idx},...
            log_std{idx}] = ...
            C.pc_signal_variance;
        
        n = sum(~isnan(signal_ipr{idx}));
        S_IPR(idx) = nanmean(signal_ipr{idx});
        S_IPR_sem(idx) = nanstd(signal_ipr{idx}) ./ sqrt(n);
        
        S_IPR_s(idx) = nanmean(signal_ipr_s{idx});
        S_IPR_s_sem(idx) = nanstd(signal_ipr_s{idx}) ./ sqrt(n);
        
        S_KL(idx) = nanmean(signal_KL{idx});
        S_KL_sem(idx) = nanstd(signal_KL{idx}) ./ sqrt(n);
        
        S_PERP(idx) = nanmean(signal_perplexity{idx});
        S_PERP_sem(idx) = nanstd(signal_perplexity{idx}) ./ sqrt(n);
        
        AREA_COS(idx) = nanmean(area_cos{idx});
        AREA_COS_sem(idx) = nanstd(area_cos{idx}) ./ sqrt(n);
        
        AREA_COS2(idx) = nanmean(area_cos2{idx});
        AREA_COS2_sem(idx) = nanstd(area_cos2{idx}) ./ sqrt(n);
        
        DMU_IPR(idx) = nanmean(dmu_ipr{idx});
        DMU_IPR_sem(idx) = nanstd(dmu_ipr{idx}) ./ sqrt(n);
        
        SD_UNBIASED(idx) = nanmean(sig_dens_unbiased{idx});
        SD_UNBIASED_sem(idx) = nanstd(sig_dens_unbiased{idx}) ./ sqrt(n);
        
        LOG_STD(idx) = nanmean(log_std{idx});
        LOG_STD_sem(idx) = nanstd(log_std{idx}) ./ sqrt(n);
        
        V_IPR(idx) = mean(spectrum_ipr{idx});
        V_IPR_sem(idx) = sem(spectrum_ipr{idx});
        progressbar(idx/numel(m_));
    end
    save(IPR_savefile, 'S_IPR', 'S_IPR_sem', 'S_IPR_s', 'S_IPR_s_sem',...
        'S_KL', 'S_KL_sem', 'S_PERP', 'S_PERP_sem', 'V_IPR', 'V_IPR_sem',...
        'AREA_COS', 'AREA_COS_sem', 'AREA_COS2', 'AREA_COS2_sem',...
        'DMU_IPR', 'DMU_IPR_sem', 'num_neurons', 'SD_UNBIASED', 'SD_UNBIASED_sem',...
        'LOG_STD', 'LOG_STD_sem');
else
    load(IPR_savefile);
end

frac = 0.5;

S_IPR_conf = 1.96*S_IPR_sem;
S_IPR_s_conf = 1.96*S_IPR_s_sem;
S_KL_conf = 1.96*S_KL_sem;
S_PERP_conf = 1.96*S_PERP_sem;
V_IPR_conf = 1.96*V_IPR_sem;
AREA_COS_conf = 1.96*AREA_COS_sem;
AREA_COS2_conf = 1.96*AREA_COS2_sem;
DMU_IPR_conf = 1.96*DMU_IPR_sem;
SD_UNBIASED_conf = 1.96*SD_UNBIASED_sem;
LOG_STD_conf = 1.96*LOG_STD_sem;

condition_S_IPR = (S_IPR_conf < frac*S_IPR);
condition_V_IPR = (V_IPR_conf < frac*V_IPR);
condition_S_KL = (S_KL_conf < frac*S_KL);
condition_A_RAT = (asymp_ratio_conf < frac*asymp_ratio);
fprintf('condS: %d failed\ncondV: %d failed\ncondKL: %d\ncondA: %d failed\n',...
    sum(~condition_S_IPR), sum(~condition_V_IPR), sum(~condition_S_KL), sum(~condition_A_RAT));

good_filter = condition_S_IPR & condition_S_KL &...
    condition_V_IPR & ...
    condition_A_RAT;

%use this to cancel the filtering of sessions, use all of them
if true
    good_filter = true | good_filter; %make it all true
end

fprintf('Only using %d out of %d sessions\n', sum(good_filter), numel(good_filter));
g_ = good_filter;

dotsize = 10;
figure;
subplot(2,2,1);
fprintf('TEST FOR SIGNAL_IPR\n');
PanelGenerator.plot_regress(S_IPR(g_), asymp_ratio(g_),...
    S_IPR_conf(g_), asymp_ratio_conf(g_), mouse_names(g_), 'k', 'dotsize', dotsize);
title(sprintf('All sessions (n=%d)', sum(g_)));
xlabel 'Signal IPR'
ylabel 'Asymp SNR / Asymp SNR shuf.'
ylim([0,Inf]);

subplot(2,2,2);
fprintf('TEST FOR VARIANCE_IPR\n');
PanelGenerator.plot_regress(V_IPR(g_), asymp_ratio(g_),...
    V_IPR_conf(g_), asymp_ratio_conf(g_), mouse_names(g_), 'k', 'dotsize', dotsize);
title(sprintf('All sessions (n=%d)', sum(g_)));
xlabel 'Variance IPR'
ylabel 'Asymp SNR / Asymp SNR shuf.'
ylim([0,Inf]);

subplot(2,2,3);
fprintf('TEST FOR SIGNAL_IPR (per mouse)\n');
PanelGenerator.plot_regress_averaged(S_IPR(g_), asymp_ratio(g_),...
    S_IPR_conf(g_), asymp_ratio_conf(g_), mouse_names(g_), 'k', 'dotsize', dotsize);
title(sprintf('Per mouse (n=%d)', numel(unique(mouse_names(g_)))));
xlabel 'Signal IPR'
ylabel 'Asymp SNR / Asymp SNR shuf.'
ylim([0,Inf]);

subplot(2,2,4);
fprintf('TEST FOR VARIANCE_IPR (per mouse)\n');
PanelGenerator.plot_regress_averaged(V_IPR(g_), asymp_ratio(g_),...
    V_IPR_conf(g_), asymp_ratio_conf(g_), mouse_names(g_), 'k', 'dotsize', dotsize);
title(sprintf('Per mouse (n=%d)', numel(unique(mouse_names(g_)))));
xlabel 'Variance IPR'
ylabel 'Asymp SNR / Asymp SNR shuf.'
ylim([0,Inf]);




figure;
subplot(2,3,1);
fprintf('TEST FOR SIGNAL_IPR_shuf\n');
PanelGenerator.plot_regress(S_IPR_s(g_), asymp_ratio(g_),...
    S_IPR_s_conf(g_), asymp_ratio_conf(g_), mouse_names(g_), 'k', 'dotsize', dotsize);
title(sprintf('All sessions (n=%d)', sum(g_)));
xlabel 'Signal IPR (shuf)'
ylabel 'Asymp SNR / Asymp SNR shuf.'
ylim([0,Inf]);

subplot(2,3,2);
fprintf('TEST FOR SIGNAL_KL\n');
PanelGenerator.plot_regress(S_KL(g_), asymp_ratio(g_),...
    S_KL_conf(g_), asymp_ratio_conf(g_), mouse_names(g_), 'k', 'dotsize', dotsize);
title(sprintf('All sessions (n=%d)', sum(g_)));
xlabel 'Signal KL'
ylabel 'Asymp SNR / Asymp SNR shuf.'
ylim([0,Inf]);

subplot(2,3,3);
fprintf('TEST FOR SIGNAL_PERP\n');
PanelGenerator.plot_regress(S_PERP(g_), asymp_ratio(g_),...
    S_PERP_conf(g_), asymp_ratio_conf(g_), mouse_names(g_), 'k', 'dotsize', dotsize);
title(sprintf('All sessions (n=%d)', sum(g_)));
xlabel 'Signal perp.'
ylabel 'Asymp SNR / Asymp SNR shuf.'
ylim([0,Inf]);

subplot(2,3,4);
fprintf('TEST FOR SIGNAL_IPR_shuf (per mouse)\n');
PanelGenerator.plot_regress_averaged(S_IPR_s(g_), asymp_ratio(g_),...
    S_IPR_s_conf(g_), asymp_ratio_conf(g_), mouse_names(g_), 'k', 'dotsize', dotsize);
title(sprintf('Per mouse (n=%d)', numel(unique(mouse_names(g_)))));
xlabel 'Signal IPR (shuf)'
ylabel 'Asymp SNR / Asymp SNR shuf.'
ylim([0,Inf]);

subplot(2,3,5);
fprintf('TEST FOR SIGNAL_KL (per mouse)\n');
PanelGenerator.plot_regress_averaged(S_KL(g_), asymp_ratio(g_),...
    S_KL_conf(g_), asymp_ratio_conf(g_), mouse_names(g_), 'k', 'dotsize', dotsize);
title(sprintf('Per mouse (n=%d)', numel(unique(mouse_names(g_)))));
xlabel 'Signal KL'
ylabel 'Asymp SNR / Asymp SNR shuf.'
ylim([0,Inf]);

subplot(2,3,6);
fprintf('TEST FOR SIGNAL_PERP (per mouse)\n');
PanelGenerator.plot_regress_averaged(S_PERP(g_), asymp_ratio(g_),...
    S_PERP_conf(g_), asymp_ratio_conf(g_), mouse_names(g_), 'k', 'dotsize', dotsize);
title(sprintf('Per mouse (n=%d)', numel(unique(mouse_names(g_)))));
xlabel 'Signal perp.'
ylabel 'Asymp SNR / Asymp SNR shuf.'
ylim([0,Inf]);

figure;
subplot(2,3,1);
fprintf('TEST FOR AREA_COS\n');
PanelGenerator.plot_regress(AREA_COS(g_), asymp_ratio(g_),...
    AREA_COS_conf(g_), asymp_ratio_conf(g_), mouse_names(g_), 'k', 'dotsize', dotsize);
title(sprintf('All sessions (n=%d)', sum(g_)));
xlabel 'Area between cos'
ylabel 'Asymp SNR / Asymp SNR shuf.'
ylim([0, Inf]);

subplot(2,3,2);
fprintf('TEST FOR AREA_COS2\n');
PanelGenerator.plot_regress(AREA_COS2(g_), asymp_ratio(g_),...
    AREA_COS2_conf(g_), asymp_ratio_conf(g_), mouse_names(g_), 'k', 'dotsize', dotsize);
title(sprintf('All sessions (n=%d)', sum(g_)));
xlabel 'Area between cos^2'
ylabel 'Asymp SNR / Asymp SNR shuf.'
ylim([0, Inf]);

subplot(2,3,3);
fprintf('TEST FOR DMU_IPR\n');
PanelGenerator.plot_regress(DMU_IPR(g_), asymp_ratio(g_),...
    DMU_IPR_conf(g_), asymp_ratio_conf(g_), mouse_names(g_), 'k', 'dotsize', dotsize);
title(sprintf('All sessions (n=%d)', sum(g_)));
xlabel 'IPR of \Delta\mu'
ylabel 'Asymp SNR / Asymp SNR shuf.'
ylim([0, Inf]);

subplot(2,3,4);
fprintf('TEST FOR AREA_COS (per mouse)\n');
PanelGenerator.plot_regress_averaged(AREA_COS(g_), asymp_ratio(g_),...
    AREA_COS_conf(g_), asymp_ratio_conf(g_), mouse_names(g_), 'k', 'dotsize', dotsize);
title(sprintf('Per mouse (n=%d)', numel(unique(mouse_names(g_)))));
xlabel 'Area between cos'
ylabel 'Asymp SNR / Asymp SNR shuf.'
ylim([0, Inf]);

subplot(2,3,5);
fprintf('TEST FOR AREA_COS2 (per mouse)\n');
PanelGenerator.plot_regress_averaged(AREA_COS2(g_), asymp_ratio(g_),...
    AREA_COS2_conf(g_), asymp_ratio_conf(g_), mouse_names(g_), 'k', 'dotsize', dotsize);
title(sprintf('Per mouse (n=%d)', numel(unique(mouse_names(g_)))));
xlabel 'Area between cos^2'
ylabel 'Asymp SNR / Asymp SNR shuf.'
ylim([0, Inf]);

subplot(2,3,6);
fprintf('TEST FOR DMU_IPR (per mouse)\n');
PanelGenerator.plot_regress_averaged(DMU_IPR(g_), asymp_ratio(g_),...
    DMU_IPR_conf(g_), asymp_ratio_conf(g_), mouse_names(g_), 'k', 'dotsize', dotsize);
title(sprintf('Per mouse (n=%d)', numel(unique(mouse_names(g_)))));
xlabel 'IPR of \Delta\mu'
ylabel 'Asymp SNR / Asymp SNR shuf.'
ylim([0, Inf]);

figure;
subplot(2,3,1);
fprintf('TEST FOR DMU_IPR\n');
PanelGenerator.plot_regress(DMU_IPR(g_), asymp_ratio(g_),...
    DMU_IPR_conf(g_), asymp_ratio_conf(g_), mouse_names(g_), 'k', 'dotsize', dotsize);
title(sprintf('All sessions (n=%d)', sum(g_)));
xlabel 'IPR of \Delta\mu'
ylabel 'Asymp SNR / Asymp SNR shuf.'
ylim([0, Inf]);

subplot(2,3,4);
fprintf('TEST FOR DMU_IPR (per mouse)\n');
PanelGenerator.plot_regress_averaged(DMU_IPR(g_), asymp_ratio(g_),...
    DMU_IPR_conf(g_), asymp_ratio_conf(g_), mouse_names(g_), 'k', 'dotsize', dotsize);
title(sprintf('Per mouse (n=%d)', numel(unique(mouse_names(g_)))));
xlabel 'IPR of \Delta\mu'
ylabel 'Asymp SNR / Asymp SNR shuf.'
ylim([0, Inf]);

subplot(2,3,2);
fprintf('TEST FOR numneurons\n');
PanelGenerator.plot_regress(num_neurons(g_), asymp_ratio(g_),...
    0*DMU_IPR_conf(g_), asymp_ratio_conf(g_), mouse_names(g_), 'k', 'dotsize', dotsize);
title(sprintf('All sessions (n=%d)', sum(g_)));
xlabel 'Number of neurons'
ylabel 'Asymp SNR / Asymp SNR shuf.'
ylim([0, Inf]);

subplot(2,3,5);
fprintf('TEST FOR numneurons (per mouse)\n');
PanelGenerator.plot_regress_averaged(num_neurons(g_), asymp_ratio(g_),...
    0*DMU_IPR_conf(g_), asymp_ratio_conf(g_), mouse_names(g_), 'k', 'dotsize', dotsize);
title(sprintf('Per mouse (n=%d)', numel(unique(mouse_names(g_)))));
xlabel 'Number of neurons'
ylabel 'Asymp SNR / Asymp SNR shuf.'
ylim([0, Inf]);


subplot(2,3,3);
fprintf('TEST FOR DMU_IPR/N\n');
PanelGenerator.plot_regress(DMU_IPR(g_)./num_neurons(g_), asymp_ratio(g_),...
    DMU_IPR_conf(g_)./num_neurons(g_), asymp_ratio_conf(g_), mouse_names(g_), 'k', 'dotsize', dotsize);
title(sprintf('All sessions (n=%d)', sum(g_)));
xlabel(sprintf('IPR of \\Delta\\mu / #neurons\n(Signal Density)'));
ylabel 'Asymp SNR / Asymp SNR shuf.'
ylim([0, Inf]);

subplot(2,3,6);
fprintf('TEST FOR DMU_IPR/N (per mouse)\n');
PanelGenerator.plot_regress_averaged(DMU_IPR(g_)./num_neurons(g_), asymp_ratio(g_),...
    DMU_IPR_conf(g_)./num_neurons(g_), asymp_ratio_conf(g_), mouse_names(g_), 'k', 'dotsize', dotsize);
title(sprintf('Per mouse (n=%d)', numel(unique(mouse_names(g_)))));
xlabel(sprintf('IPR of \\Delta\\mu / #neurons\n(Signal Density)'));
ylabel 'Asymp SNR / Asymp SNR shuf.'
ylim([0, Inf]);

figure; %for pablo...
subplot(2,1,1);
fprintf('TEST FOR AREA_COS/N to the inverse ratio\n');
PanelGenerator.plot_regress(AREA_COS(g_)./num_neurons(g_), 1./asymp_ratio(g_),...
    AREA_COS_conf(g_)./num_neurons(g_), asymp_ratio_conf(g_)./asymp_ratio(g_).^2, mouse_names(g_), 'k', 'dotsize', dotsize);
title(sprintf('All sessions (n=%d)', sum(g_)));
xlabel 'Area between cos / #neurons'
ylabel 'Asymp SNR shuf./ Asymp SNR'
ylim([0, Inf]);

subplot(2,1,2);
fprintf('TEST FOR AREA_COS/N to the inverse ratio (per mouse)\n');
PanelGenerator.plot_regress_averaged(AREA_COS(g_)./num_neurons(g_), 1./asymp_ratio(g_),...
    AREA_COS_conf(g_)./num_neurons(g_), asymp_ratio_conf(g_)./asymp_ratio(g_).^2, mouse_names(g_), 'k', 'dotsize', dotsize);
title(sprintf('Per mouse (n=%d)', numel(unique(mouse_names(g_)))));
xlabel 'Area between cos / #neurons'
ylabel 'Asymp SNR shuf./ Asymp SNR'
ylim([0, Inf]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
subplot(2,5,1);
fprintf('TEST FOR AREA_COS\n');
PanelGenerator.plot_regress(AREA_COS(g_), asymp_ratio(g_),...
    AREA_COS_conf(g_), asymp_ratio_conf(g_), mouse_names(g_), 'k', 'dotsize', dotsize);
title(sprintf('All sessions (n=%d)', sum(g_)));
xlabel 'Area between cos'
ylabel 'Asymp SNR / Asymp SNR shuf.'
ylim([0, 1]);

subplot(2,5,2);
fprintf('TEST FOR AREA_COS2\n');
PanelGenerator.plot_regress(AREA_COS2(g_), asymp_ratio(g_),...
    AREA_COS2_conf(g_), asymp_ratio_conf(g_), mouse_names(g_), 'k', 'dotsize', dotsize);
title(sprintf('All sessions (n=%d)', sum(g_)));
xlabel 'Area between cos^2'
ylabel 'Asymp SNR / Asymp SNR shuf.'
ylim([0, 1]);


subplot(2,5,6);
fprintf('TEST FOR AREA_COS (per mouse)\n');
PanelGenerator.plot_regress_averaged(AREA_COS(g_), asymp_ratio(g_),...
    AREA_COS_conf(g_), asymp_ratio_conf(g_), mouse_names(g_), 'k', 'dotsize', dotsize);
title(sprintf('Per mouse (n=%d)', numel(unique(mouse_names(g_)))));
xlabel 'Area between cos'
ylabel 'Asymp SNR / Asymp SNR shuf.'
ylim([0, 1]);

subplot(2,5,7);
fprintf('TEST FOR AREA_COS2 (per mouse)\n');
PanelGenerator.plot_regress_averaged(AREA_COS2(g_), asymp_ratio(g_),...
    AREA_COS2_conf(g_), asymp_ratio_conf(g_), mouse_names(g_), 'k', 'dotsize', dotsize);
title(sprintf('Per mouse (n=%d)', numel(unique(mouse_names(g_)))));
xlabel 'Area between cos^2'
ylabel 'Asymp SNR / Asymp SNR shuf.'
ylim([0, 1]);

subplot(2,5,3);
fprintf('TEST FOR numneurons\n');
PanelGenerator.plot_regress(num_neurons(g_), asymp_ratio(g_),...
    0*DMU_IPR_conf(g_), asymp_ratio_conf(g_), mouse_names(g_), 'k', 'dotsize', dotsize);
title(sprintf('All sessions (n=%d)', sum(g_)));
xlabel 'Number of neurons'
ylabel 'Asymp SNR / Asymp SNR shuf.'
ylim([0, 1]);

subplot(2,5,8);
fprintf('TEST FOR numneurons (per mouse)\n');
PanelGenerator.plot_regress_averaged(num_neurons(g_), asymp_ratio(g_),...
    0*DMU_IPR_conf(g_), asymp_ratio_conf(g_), mouse_names(g_), 'k', 'dotsize', dotsize);
title(sprintf('Per mouse (n=%d)', numel(unique(mouse_names(g_)))));
xlabel 'Number of neurons'
ylabel 'Asymp SNR / Asymp SNR shuf.'
ylim([0, 1]);

subplot(2,5,4);
fprintf('TEST FOR DMU_IPR/N\n');
PanelGenerator.plot_regress(DMU_IPR(g_)./num_neurons(g_), asymp_ratio(g_),...
    DMU_IPR_conf(g_)./num_neurons(g_), asymp_ratio_conf(g_), mouse_names(g_), 'k', 'dotsize', dotsize);
title(sprintf('All sessions (n=%d)', sum(g_)));
xlabel(sprintf('IPR of \\Delta\\mu / #neurons\n(Signal Density)'));
ylabel 'Asymp SNR / Asymp SNR shuf.'
ylim([0, Inf]);

subplot(2,5,9);
fprintf('TEST FOR DMU_IPR/N (per mouse)\n');
PanelGenerator.plot_regress_averaged(DMU_IPR(g_)./num_neurons(g_), asymp_ratio(g_),...
    DMU_IPR_conf(g_)./num_neurons(g_), asymp_ratio_conf(g_), mouse_names(g_), 'k', 'dotsize', dotsize);
title(sprintf('Per mouse (n=%d)', numel(unique(mouse_names(g_)))));
xlabel(sprintf('IPR of \\Delta\\mu / #neurons\n(Signal Density)'));
ylabel 'Asymp SNR / Asymp SNR shuf.'
ylim([0, Inf]);

%for pablo...
subplot(2,5,5);
fprintf('TEST FOR AREA_COS/N to the inverse ratio\n');
PanelGenerator.plot_regress(AREA_COS(g_)./num_neurons(g_), 1./asymp_ratio(g_),...
    AREA_COS_conf(g_)./num_neurons(g_), asymp_ratio_conf(g_)./asymp_ratio(g_).^2, mouse_names(g_), 'k', 'dotsize', dotsize);
title(sprintf('All sessions (n=%d)', sum(g_)));
xlabel 'Area between cos / #neurons'
ylabel 'Asymp SNR shuf./ Asymp SNR'
ylim([0, Inf]);

subplot(2,5,10);
fprintf('TEST FOR AREA_COS/N to the inverse ratio (per mouse)\n');
PanelGenerator.plot_regress_averaged(AREA_COS(g_)./num_neurons(g_), 1./asymp_ratio(g_),...
    AREA_COS_conf(g_)./num_neurons(g_), asymp_ratio_conf(g_)./asymp_ratio(g_).^2, mouse_names(g_), 'k', 'dotsize', dotsize);
title(sprintf('Per mouse (n=%d)', numel(unique(mouse_names(g_)))));
xlabel 'Area between cos / #neurons'
ylabel 'Asymp SNR shuf./ Asymp SNR'
ylim([0, Inf]);


figure;
subplot(2,3,1);
fprintf('TEST FOR DMU_IPR/N\n');
PanelGenerator.plot_regress(DMU_IPR(g_)./num_neurons(g_), asymp_ratio(g_),...
    DMU_IPR_conf(g_)./num_neurons(g_), asymp_ratio_conf(g_), mouse_names(g_), 'k', 'dotsize', dotsize);
title(sprintf('All sessions (n=%d)', sum(g_)));
xlabel(sprintf('IPR of \\Delta\\mu / #neurons\n(Signal Density)'));
ylabel 'Asymp SNR / Asymp SNR shuf.'
ylim([0, Inf]);

subplot(2,3,4);
fprintf('TEST FOR DMU_IPR/N (per mouse)\n');
PanelGenerator.plot_regress_averaged(DMU_IPR(g_)./num_neurons(g_), asymp_ratio(g_),...
    DMU_IPR_conf(g_)./num_neurons(g_), asymp_ratio_conf(g_), mouse_names(g_), 'k', 'dotsize', dotsize);
title(sprintf('Per mouse (n=%d)', numel(unique(mouse_names(g_)))));
xlabel(sprintf('IPR of \\Delta\\mu / #neurons\n(Signal Density)'));
ylabel 'Asymp SNR / Asymp SNR shuf.'
ylim([0, Inf]);

subplot(2,3,2);
fprintf('TEST FOR SD_UNBIASED\n');
PanelGenerator.plot_regress(SD_UNBIASED(g_), asymp_ratio(g_),...
    SD_UNBIASED_conf(g_), asymp_ratio_conf(g_), mouse_names(g_), 'k', 'dotsize', dotsize);
title(sprintf('All sessions (n=%d)', sum(g_)));
xlabel(sprintf('Unbiased signal density'));
ylabel 'Asymp SNR / Asymp SNR shuf.'
ylim([0, Inf]);

subplot(2,3,5);
fprintf('TEST FOR SD_UNBIASED (per mouse)\n');
PanelGenerator.plot_regress_averaged(SD_UNBIASED(g_), asymp_ratio(g_),...
    SD_UNBIASED_conf(g_), asymp_ratio_conf(g_), mouse_names(g_), 'k', 'dotsize', dotsize);
title(sprintf('Per mouse (n=%d)', numel(unique(mouse_names(g_)))));
xlabel(sprintf('Unbiased signal density'));
ylabel 'Asymp SNR / Asymp SNR shuf.'
ylim([0, Inf]);

subplot(2,3,3);
fprintf('TEST FOR LOG_STD\n');
PanelGenerator.plot_regress(LOG_STD(g_), asymp_ratio(g_),...
    LOG_STD_conf(g_), asymp_ratio_conf(g_), mouse_names(g_), 'k', 'dotsize', dotsize);
title(sprintf('All sessions (n=%d)', sum(g_)));
xlabel(sprintf('Log stdev.'));
ylabel 'Asymp SNR / Asymp SNR shuf.'
ylim([0, Inf]);

subplot(2,3,6);
fprintf('TEST FOR LOG_STD (per mouse)\n');
PanelGenerator.plot_regress_averaged(LOG_STD(g_), asymp_ratio(g_),...
    LOG_STD_conf(g_), asymp_ratio_conf(g_), mouse_names(g_), 'k', 'dotsize', dotsize);
title(sprintf('Per mouse (n=%d)', numel(unique(mouse_names(g_)))));
xlabel(sprintf('Log stdev.'));
ylabel 'Asymp SNR / Asymp SNR shuf.'
ylim([0, Inf]);
end


function [quotient, quotient_uncertainty] = uncertain_divide(x, xc, y, yc)
quotient = x./y;
quotient_uncertainty = abs(x./y).*sqrt((xc./x).^2 + (yc./y).^2);
end