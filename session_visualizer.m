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
    end
    
    methods
        function sv = session_visualizer(o)
            sv.o = o;
            sv.fig = figure;
            
            
            subplot(2,2,1);
            sv.y_plot = plot(sv.o.data.y.raw.full); hold on;
            sv.y_plot.ButtonDownFcn = @(src, event) sv.trial_finder(event.IntersectionPoint);
            sv.y_plot_trial = plot(sv.o.data.y.raw.full);
            
            
            subplot(2,2,3);
            [~, ~, sv.act,~,~,~,~,sv.stats] = plsregress(o.data.X.fast, o.data.y.scaled, 2);
            sv.xmean = mean(sv.o.data.X.fast);
            sv.pls_plot = scatter(sv.act(:,1), sv.act(:,2),...
                20, sv.o.data.y.raw.fast, 'filled', 'MarkerFaceAlpha', 0.05); hold on;
            sv.pls_plot.ButtonDownFcn = @(src, event) sv.scatter_trial_finder(event.IntersectionPoint);
            sv.pls_plot_trial = scatter(sv.act(:,1), sv.act(:,2), 50, sv.o.data.y.raw.fast, 'filled');
            sv.trial_index = 100;
            
            
            subplot(2,2,[2 4]);
            sv.trial_level_pca;
            sv.trial_pca.ButtonDownFcn = @(src, event) sv.trial_pca_picker(event.IntersectionPoint);
            sv.trial_pca_trial = scatter(sv.tpca(:,1), sv.tpca(:,2), 60, 'r');
            
            sv.refresh;
        end
        
        function refresh(sv)
            sv.trial_pca_trial.XData = sv.tpca(sv.trial_index, 1);
            sv.trial_pca_trial.YData = sv.tpca(sv.trial_index, 2);
            
            sel = sv.trial_start(sv.trial_index):sv.trial_end(sv.trial_index);
            sv.y_plot_trial.XData = sel;
            sv.y_plot_trial.YData = sv.o.data.y.raw.full(sel);
            
            sel = sv.trial_start(sv.trial_index):sv.trial_end(sv.trial_index);
            dat = (sv.o.data.X.full(sel,:) - sv.xmean) * sv.stats.W;
            sv.pls_plot_trial.XData = dat(:,1);
            sv.pls_plot_trial.YData = dat(:,2);
            sv.pls_plot_trial.CData = sv.o.data.y.raw.full(sel);
        end
        
        
        function trial_pca_picker(sv, ip)
            tpca = [sv.trial_pca.XData.', sv.trial_pca.YData.'];
            d = sum((ip(1:2) - tpca).^2,2);
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
            
            sv.trial_pca = scatter(score_(:,1), score_(:,2), 20, trial_end - trial_start, 'filled');
            sv.tpca = score_;
            hold on; trial_origin = -mean(trial_feature_flat)*coeff_;
            scatter(trial_origin(1), trial_origin(2), 40, 'r', 'filled');
            title('Trial-level PCA, color is length');
            
        end
        
        
        function trial_finder(sv, ip)
            f_ind = round(ip(1));
            tr_ids = cumsum(diff(sv.o.data.mask.fast) == 1);
            sv.trial_index = tr_ids(f_ind);
            sv.refresh;
        end
        
        function scatter_trial_finder(sv, ip)
            d = sum((ip(1:2) - sv.act).^2,2);
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