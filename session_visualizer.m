classdef session_visualizer < handle
    properties
        o;
        act;
        y_plot;
        y_plot_trial;
        pls_plot;
        pls_plot_trial;
        fig;
        trial_start;
        trial_end;
        stats;
        trial_index;
        xmean;
        trial_pca;
        trial_pca_trial;
        tpca;
        trial_raster;
        fw_ord;
        bw_ord;
        ks_plot;
        preds;
        pred_plot;
    end
    
    methods
        function sv = session_visualizer(o)
            sv.o = o;
            sv.fig = figure;
            
            
            subplot(4,4,[1 2]);
            sv.y_plot = plot(sv.o.data.y.raw.full, 'DisplayName', 'Position trace'); hold on;
            sv.y_plot.ButtonDownFcn = @(src, event) sv.trial_finder(event.IntersectionPoint);
            sv.y_plot_trial = plot(sv.o.data.y.raw.full, 'DisplayName', 'Selected trial');
            legend;
            xlabel('Frames (20Hz)');
            ylabel('Position (px)');
            title('Measured mouse position');
            
            subplot(4,4,[5 6]);
            sv.ks_plot = plot(sv.o.data.y.ks, '-or', 'DisplayName', 'True place bin');
            hold on;
            sv.preds = mean([o.res.unshuf.te_pred{end,:}],2);
            sv.pred_plot = plot(sv.preds, '-ob', 'DisplayName', 'Prediction');
            ylim([1 2*sv.o.opt.n_bins]);
            legend;
            xlabel('Frames (20Hz)');
            ylabel('Position bin (f/w odd, b/w even)');
            title('Place decoding on selected trial');
            
            subplot(4,4,[9 10 13 14]);
            [~, ~, sv.act,~,~,~,~,sv.stats] = plsregress(o.data.X.fast, zscore([o.data.y.scaled, o.data.y.direction]), 2);
            sv.xmean = mean(sv.o.data.X.fast);
            sv.pls_plot = scatter(sv.act(:,1), sv.act(:,2),...
                20, sv.o.data.y.scaled, 'filled', 'MarkerFaceAlpha', 0.05, 'DisplayName', 'All trials'); hold on;
            sv.pls_plot.ButtonDownFcn = @(src, event) sv.scatter_trial_finder(event.IntersectionPoint);
            sv.pls_plot_trial = scatter(sv.act(:,1), sv.act(:,2), 50, sv.o.data.y.scaled, 'filled', 'DisplayName', 'Selected trial');
            pls_origin = -sv.xmean * sv.stats.W;
            scatter(pls_origin(1), pls_origin(2), 100, 'r', 'DisplayName', 'Origin');
            legend;
            xlabel('PLS1'); ylabel('PLS2');
            title('2D PLS regression on track position');
            colorbar;
            sv.trial_index = 1;
            
            
            subplot(4,4,[3 7 11 15]);
            sv.trial_level_pca;
            sv.trial_pca.ButtonDownFcn = @(src, event) sv.trial_pca_picker(event.IntersectionPoint);
            sv.trial_pca_trial = scatter(sv.tpca(:,2), sv.tpca(:,1), 60, 'r');
            xlabel('PC2'); ylabel('PC1'); colorbar;
            
            subplot(4,4,[4 8 12 16]);
            sv.trial_raster = imagesc(sv.o.data.X.fast);
            xlabel('Neurons, position ordered');
            ylabel('Frames (20Hz)');
            
            q_ = sv.o.bin_data(false, false);
            mean_bin_X = q_.bin_X.mean;
            fw_bindist = mean_bin_X(1:2:end,:) ./ sum(mean_bin_X(1:2:end,:),1);
            bw_bindist = mean_bin_X(2:2:end,:) ./ sum(mean_bin_X(2:2:end,:),1);
            fw_mb = (1:o.opt.n_bins)*fw_bindist;
            bw_mb = (1:o.opt.n_bins)*bw_bindist;
            [~, sv.fw_ord] = sort(fw_mb);
            [~, sv.bw_ord] = sort(bw_mb);
            
            colormap jet;
            sv.refresh;
        end
        
        function refresh(sv)
            sel = sv.trial_start(sv.trial_index):sv.trial_end(sv.trial_index);
            
            trans = cumsum(sv.o.data.mask.fast);
            
            sv.ks_plot.XData = 1:numel(sel);
            sv.ks_plot.YData = sv.o.data.y.ks(trans(sel));
            sv.pred_plot.XData = 1:numel(sel);
            sv.pred_plot.YData = sv.preds(trans(sel));
            
            sv.trial_pca_trial.XData = sv.tpca(sv.trial_index, 2);
            sv.trial_pca_trial.YData = sv.tpca(sv.trial_index, 1);
            
            
            sv.y_plot_trial.XData = sel;
            sv.y_plot_trial.YData = sv.o.data.y.raw.full(sel);
            
            dat = (sv.o.data.X.full(sel,:) - sv.xmean) * sv.stats.W;
            sv.pls_plot_trial.XData = dat(:,1);
            sv.pls_plot_trial.YData = dat(:,2);
            sv.pls_plot_trial.CData = sv.o.data.y.scaled(trans(sel));
            sv.pls_plot_trial.DisplayName = sprintf('Trial %d', sv.trial_index);
            
            if round(mean(sv.o.data.y.direction(trans(sel)))) == 1
                ord_ = sv.fw_ord;
                mes_ = 'Neural activity, neurons ordered by fw direction selectivity';
            else
                ord_ = sv.bw_ord;
                mes_ = 'Neural activity, neurons ordered by bw direction selectivity';
            end
            raster = sv.o.data.X.full(sel,ord_);
            subplot(4,4,[4 8 12 16]);
            sv.trial_raster = imagesc(flipud(raster.'));
            ylabel('Neurons, position ordered');
            xlabel('Frames (20Hz)');
            title(mes_);
            colorbar;
        end
        
        
        function trial_pca_picker(sv, ip)
            d = sum((ip([2 1]) - sv.tpca).^2,2);
            [~, sv.trial_index] = min(d);
            
            sv.refresh;
        end
        
        function trial_level_pca(sv)
            o = sv.o;
            msk = [0; o.data.mask.fast];
            sv.trial_start = find(diff(msk) == 1);
            sv.trial_end = find(diff(msk) == -1);
            trial_start = sv.trial_start;
            trial_end = sv.trial_end;
            
            i_convert = cumsum(o.data.mask.fast);
            num_trials = numel(trial_start);
            
            
            trial_feature = zeros(num_trials, o.opt.n_bins, o.data.total_neurons);
            trial_dir = zeros(1, num_trials);
            for t_i = 1:num_trials
                sel = i_convert(trial_start(t_i):trial_end(t_i));
                sel_save{t_i} = sel;
                trial_X = o.data.X.fast(sel,:);
                trial_ks = o.data.y.ks(sel);
                trial_ks(mod(trial_ks,2)==1) = trial_ks(mod(trial_ks,2)==1) + 1;
                trial_ks = trial_ks/2;
                for b = 1:o.opt.n_bins
                    trial_feature(t_i, b, :) = mean(trial_X(trial_ks==b,:));
                end
                trial_dir(t_i) = mean(o.data.y.direction(sel));
            end
            
            trial_feature(isnan(trial_feature)) = 0;
            trial_feature_flat = reshape(trial_feature, [], o.opt.n_bins*o.data.total_neurons);
            [coeff_, score_, ~] = pca(trial_feature_flat, 'NumComponents', 2);
            
            sv.trial_pca = scatter(score_(:,2), score_(:,1), 20, trial_end - trial_start, 'filled');
            sv.tpca = score_;
            hold on; trial_origin = -mean(trial_feature_flat)*coeff_;
            scatter(trial_origin(2), trial_origin(1), 40, 'r', 'filled');
            title('Trial-level PCA, color is length');
            
        end
        
        
        function trial_finder(sv, ip)
            f_ind = round(ip(1));
            tr_ids = cumsum(diff(sv.o.data.mask.fast) == 1);
            sv.trial_index = tr_ids(f_ind);
            sv.refresh;
        end
        
        function scatter_trial_finder(sv, ip)
            d = sum((ip([1 2]) - sv.act).^2,2);
            [~, ind] = min(d);
            qw = cumsum(sv.o.data.mask.fast);
            f_ind = find(qw == ind, 1);
            
            tr_ids = cumsum(diff(sv.o.data.mask.fast) == 1);
            sv.trial_index = tr_ids(f_ind);
            sv.refresh;
        end
    end
    
end
%
% f = figure;
% subplot(2,1,1);
% p1 = plot(o.data.y.raw.full);
%
% subplot(2,1,2);
% [XL, ~, act,~,~,~,~,stats] = plsregress(o.data.X.fast, o.data.y.scaled, 2);
% p2 = scatter(act(:,1), act(:,2), 4, o.data.y.scaled);
%
%
% trial_start = find(diff(o.data.mask.fast) == 1);
% trial_end = find(diff(o.data.mask.fast) == -1);
%
% tr_ix = 100;
% refresh(tr_ix, act, stats, o, trial_start, trial_end);
%
%
% function search_for_point(ip, act, stats, o, trial_start, trial_end)
% d = sum((ip(1:2) - act).^2,2);
% [~, ind] = min(d);
% qw = cumsum(o.data.mask.fast);
% f_ind = find(qw == ind, 1);
%
% tr_ids = cumsum(diff(o.data.mask.fast) == 1);
% tr_ix = tr_ids(f_ind);
%
% refresh(tr_ix, act, stats, o, trial_start, trial_end);
% end
%
% function trial_finder(ip, act, stats, o, trial_start, trial_end)
% f_ind = round(ip(1));
%
% tr_ids = cumsum(diff(o.data.mask.fast) == 1);
% tr_ix = tr_ids(f_ind);
% refresh(tr_ix, act, stats, o, trial_start, trial_end);
% end
%
% function refresh(tr_ix, act, stats, o, trial_start, trial_end)
% tr_range = trial_start(tr_ix):trial_end(tr_ix);
% subplot(2,1,1);
% xl_ = xlim;
% pp = plot(o.data.y.raw.full);
% pp.ButtonDownFcn = @(src, event) trial_finder(event.IntersectionPoint, act, stats, o, trial_start, trial_end);
% hold on;
% plot(tr_range, o.data.y.raw.full(tr_range)); hold off;
% xlim(xl_);
% subplot(2,1,2);
% sp = scatter(act(:,1), act(:,2), 20, o.data.y.raw.fast,...
%     'filled', 'MarkerFaceAlpha', 0.05);
% sp.ButtonDownFcn = @(src, event) search_for_point(event.IntersectionPoint, act, stats, o, trial_start, trial_end);
% hold on;
% xmean = mean(o.data.X.fast);
% dat = (o.data.X.full(tr_range,:) - xmean) * stats.W;
% scatter(dat(:,1), dat(:,2), 50, o.data.y.raw.full(tr_range), 'filled');
% hold off;
% end