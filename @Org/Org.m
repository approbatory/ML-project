classdef Org < handle

    properties
        vars;
        derived;
        storage_file = 'default_store.mat'
        total_sessions;
        
        mouse;
        
        sess_prop = [];
        sess_prop_conf = [];
    end
    
    methods
        function o = Org
            sm = SessManager;
            o.total_sessions = sm.num_usable;
            if ~exist(o.storage_file, 'file')
                o.vars = struct;
            else
                L = load(o.storage_file);
                o.vars = L.vars_;
                o.mouse = L.mouse_;
                o.derived = L.derived_;
            end
        end
        
        function init(o)
            o.load_definitions;
            o.session_properties_dataset;
        end
        
        function save_me(o)
            if ~exist(o.storage_file, 'file')
                vars_ = o.vars;
                mouse_ = o.mouse;
                save('-v7.3', o.storage_file, 'vars_', 'mouse_');
            end
            m = matfile(o.storage_file, 'Writable', true);
            m.derived_ = o.derived;
        end
        
        res = fetch(o, varname)
        load_from_clouds(o)
        load_definitions(o)
        
        res = sess_by_bins(o, varname, sess_idx)
        res = sess_med_bins(o, varname, sess_idx)
        
        res_list = mouse_all_sess(o, varname, mouse_name)
        [res, res_sem] = mouse_by_bins(o, varname, mouse_name)
        [res, res_sem] = mouse_med_bins(o, varname, mouse_name)
        
        [res, res_sem] = all_by_bins(o, varname, restrict)
        [res, res_sem] = all_med_bins(o, varname, restrict)
        
        [res, res_sem] = per_sess(o, varname)
        
        function make_derived(o, varname, varlist, func, saveit)
            o.derived.(varname).v = varlist;
            o.derived.(varname).f = func;
            
            if exist('saveit', 'var') && saveit
                o.vars.(varname) = o.fetch(varname);
            end
        end
        
        function make_var_per_sess(o, varname, varlist, func)
            for i = 1:o.total_sessions
                for j = 1:numel(varlist)
                    my_var = varlist{j};
                    %invar{j} = o.vars.(my_var){i};
                    t_ = o.fetch(my_var);
                    invar{j} = t_{i};
                end
                outvar{i} = func(invar{:});
            end
            o.vars.(varname) = outvar;
        end
        
        [mat,sem,varnames,mat_agg,sem_agg] = predictor_matrix(o)
        [mean_table, sem_table] = session_properties_dataset(o)
        
        [T, asymp_ratio, mat, varnames] = correspondence(o, use_n50, plot_all)
        
        com_dist_vs_corr(o)
        
        function adjr2 = correlogram_old(o, var1, var2, aggregate)
            if iscell(var1)
                v1 = var1{1}; v1sem = var1{2};
            else
                [v1, v1sem] = o.per_sess(var1);
            end
            if iscell(var2)
                v2 = var2{1}; v2sem = var2{2};
            else
                [v2, v2sem] = o.per_sess(var2); 
            end
            if aggregate
                func = @PanelGenerator.plot_regress_averaged;
            else
                func = @PanelGenerator.plot_regress;
            end
            adjr2 = func(v1(:), v2(:), v1sem(:)*1.96, v2sem(:)*1.96, o.mouse, 'r', 'dotsize', 10);
            if ~iscell(var1)
                xlabel(esc(var1));
            end
            if ~iscell(var2)
                ylabel(esc(var2));
            end
        end
        
        function T = corr_table(o, var_name, highqual, conf_filt)
            if ~exist('conf_filt', 'var')
                conf_filt = false;
            end
            if ~exist('highqual', 'var')
                highqual = false;
            end
            %low_bounds = Org.default_low_bounds;
            %filt = o.sess_prop.num_neurons >= low_bounds(1) &...
            %    o.sess_prop.num_trials >= low_bounds(2);
            if highqual
                filt = SessManager.highqual_filt_from_usable;
            else
                filt = true(1,o.total_sessions);
            end
            
            [~, ~, pred_names, ~] = o.predictor_matrix;
            var_target_ = o.sess_prop.(var_name);
            var_target_conf_ = o.sess_prop_conf.(var_name);
            
            if conf_filt
                filt = filt & (var_target_' > var_target_conf_');
            end
            var_target = var_target_(filt);
            var_target_conf = var_target_conf_(filt);
            
            my_mice = o.mouse;
            my_mice = my_mice(filt);
            
            [var_target_agg, var_target_conf_agg] = Org.agg(my_mice,...
                var_target, var_target_conf);
            
            for i = 1:numel(pred_names)
                pred_var_name = pred_names{i};
                var_pred = o.sess_prop.(pred_var_name);
                var_pred = var_pred(filt);
                var_pred_conf = o.sess_prop_conf.(pred_var_name);
                var_pred_conf = var_pred_conf(filt);
                [pearson(i,1), pearson_p(i,1),...
                    spearman(i,1), spearman_p(i,1),...
                    kendall(i,1), kendall_p(i,1),...
                    adjr2(i,1)] = Org.corr_check(var_pred, var_target);
                [var_pred_agg, var_pred_conf_agg] = Org.agg(my_mice,...
                    var_pred, var_pred_conf);
                [pearson_agg(i,1), pearson_p_agg(i,1),...
                    spearman_agg(i,1), spearman_p_agg(i,1),...
                    kendall_agg(i,1), kendall_p_agg(i,1),...
                    adjr2_agg(i,1)] = Org.corr_check(var_pred_agg, var_target_agg);
            end
            T = table(pearson, pearson_p, spearman, spearman_p,...
                kendall, kendall_p, adjr2,...
                pearson_agg, pearson_p_agg, spearman_agg, spearman_p_agg,...
                kendall_agg, kendall_p_agg, adjr2_agg, 'RowNames', pred_names);
        end
        function metrics = correlogram(o, var1, var2, aggregate, highqual, conf_filt)
            if ~exist('conf_filt', 'var')
                conf_filt = false;
            end
            %if ~exist('low_bounds', 'var')
            %    low_bounds = Org.default_low_bounds;
            %end
            v1 = o.sess_prop.(var1);
            v1_conf = o.sess_prop_conf.(var1);
            
            v2 = o.sess_prop.(var2);
            v2_conf = o.sess_prop_conf.(var2);
            
            my_mice = o.mouse;
            
            %filt = o.sess_prop.N50 > 2*o.sess_prop_conf.N50;
            %filt = filt & (v1 > 2*v1_conf) & (v2 > 2*v2_conf); %<>
            %filt = o.sess_prop.num_neurons >= low_bounds(1) &...
            %    o.sess_prop.num_trials >= low_bounds(2);
            if ~exist('highqual', 'var')
                highqual = false;
            end

            if highqual
                filt = SessManager.highqual_filt_from_usable;
            else
                filt = true(1,o.total_sessions);
            end
            if conf_filt
                filt = filt & (v1' > v1_conf') & (v2' > v2_conf');
            end
            if strcmp(var2, 'asymp_ratio')
                filt = filt & (v2' < 7);
                fprintf('Special for asymp_ratio, removing outliers > 7\n');
            end
            if any(~filt)
                fprintf('Using only %d out of %d sessions\n',...
                    sum(filt), numel(filt));
                %excluded = find(~filt);
                %fprintf('Throwing out these sess indices:\n');
                %disp(excluded);
            end
            v1 = v1(filt);
            v1_conf = v1_conf(filt);
            v2 = v2(filt);
            v2_conf = v2_conf(filt);
            my_mice = my_mice(filt);
            
            if aggregate
                func = @PanelGenerator.plot_regress_averaged;
                [v1_agg, v1_conf_agg] = Org.agg(my_mice, v1, v1_conf);
                [v2_agg, v2_conf_agg] = Org.agg(my_mice, v2, v2_conf);
                [p,pp,s,sp,k,kp,adjr2] = Org.corr_check(v1_agg, v2_agg);
                fprintf('Mouse-aggregated correlations %s vs. %s: adj. R^2 = %.3f\n', var1, var2, adjr2);
                fprintf('Pearson: %.3f, p = %e, %s\n', p, pp, Utils.pstar(pp));
                fprintf('Spearman: %.3f, p = %e, %s\n', s, sp, Utils.pstar(sp));
                fprintf('Kendall: %.3f, p = %e, %s\n', k, kp, Utils.pstar(kp));
            else
                func = @PanelGenerator.plot_regress;
                [p,pp,s,sp,k,kp,adjr2] = Org.corr_check(v1, v2);
                fprintf('Sessionwise correlations %s vs. %s: adj. R^2 = %.3f\n', var1, var2, adjr2);
                fprintf('Pearson: %.3f, p = %e, %s\n', p, pp, Utils.pstar(pp));
                fprintf('Spearman: %.3f, p = %e, %s\n', s, sp, Utils.pstar(sp));
                fprintf('Kendall: %.3f, p = %e, %s\n', k, kp, Utils.pstar(kp));
            end
            adjr2_ = func(v1(:), v2(:), v1_conf(:), v2_conf(:), my_mice, 'k', 'dotsize', 10);
            assert(abs(adjr2 - adjr2_) < 0.01, 'Mismatch in adj. R^2 values');
            xlabel(esc(var1));
            ylabel(esc(var2));
            
            metrics.p = p; metrics.pp = pp; metrics.s = s; metrics.sp = sp;
            metrics.k = k; metrics.kp = kp; metrics.adjr2 = adjr2;
            metrics.aggregated = aggregate;
        end
        
        function bns(o, var1, var2, aggregate, logscale, highqual)
            if ~exist('highqual', 'var')
                highqual = false;
            end
            [v1, v1_conf, mice] = o.fetch_filt(var1, highqual);
            [v2, v2_conf, ~] = o.fetch_filt(var2, highqual);
            
            Utils.bns_groupings(v1', v2', v1_conf', v2_conf',...
                mice, aggregate, {'Real data', 'Shuffled'}, logscale);
        end
        
        function boxplot(o, var1, var2, logscale, highqual)
            [v1, v1_conf, mice] = o.fetch_filt(var1, highqual);
            [v2, v2_conf, ~] = o.fetch_filt(var2, highqual);
            
            hold on;
            bf_ = @(x, c) boxplot(x, Utils.cf_(@(x)x(end-1:end), mice), 'Colors', c, 'BoxStyle', 'outline', 'Symbol', '');
            bf_(v1, 'b'); bf_(v2, 'r');
            xtickangle(45);
            if logscale
                set(gca, 'YScale', 'log');
            end
            x_coords = get(gca, 'XTick');
            mouse_id = get(gca, 'XTickLabel');
            mouse_id = Utils.cf_(@(x)['Mouse20' x], mouse_id);
            mouse_colors = DecodeTensor.mcolor(mouse_id);
            for i = 1:numel(x_coords)
                x_val = x_coords(i);
                color = mouse_colors{i};
                hold on;
                yt_ = get(gca, 'YTick');
                %if logscale
                %    y_bottom = exp(log(yl_(1))+0.1*log(yl_(2)));
                %else
                %    y_bottom = yl_(1)+0.1*yl_(2);
                %end
                low_vals = [v1(:)-v1_conf(:); v2(:)-v2_conf(:)];
                low_vals = low_vals(low_vals > 0);
                y_bottom = min(low_vals);
                scatter(x_val, y_bottom, 8, color, 'filled');
            end
        end
        
        function [v, v_conf, mice] = fetch_filt(o, varname, highqual)
            v = o.sess_prop.(varname);
            v_conf = o.sess_prop_conf.(varname);
            
            %low_bounds = Org.default_low_bounds;
            %filt = o.sess_prop.num_neurons >= low_bounds(1) &...
            %    o.sess_prop.num_trials >= low_bounds(2);
            if ~exist('highqual', 'var')
                highqual = false;
            end
            
            if highqual
                filt = SessManager.highqual_filt_from_usable;
            else
                filt = true(1,o.total_sessions);
            end
            v = v(filt);
            v_conf = v_conf(filt);
            
            mice = o.mouse;
            mice = mice(filt);
        end
        
        function filt = default_filt(o, highqual)
            %low_bounds = Org.default_low_bounds;
            %filt = o.sess_prop.num_neurons >= low_bounds(1) &...
            %    o.sess_prop.num_trials >= low_bounds(2);
            if ~exist('highqual', 'var')
                highqual = false;
            end
            
            if highqual
                filt = SessManager.highqual_filt_from_usable;
            else
                filt = true(1,o.total_sessions);
            end
        end
        
        eigen_snr_crossval_aggregation(org, just_snr);
        area_between_cos2(o, MAX_DIM, lean);
        sigdens_plots(o);
        
        change_over_days(o, varname, highqual);
    end
    
    methods(Static)
        
        
        function [r, rc] = agg(mouse, x, xc)
            uniq_mice = unique(mouse);
            n = numel(uniq_mice);
            r = zeros(n, size(x,2));
            rc = zeros(n, size(xc,2));
            for i = 1:n
                filt = strcmp(mouse, uniq_mice{i});
                r(i,:) = mean(x(filt,:),1);
                rc(i,:) = sqrt(mean(xc(filt,:).^2,1));
            end
        end
        
        function [p,pp,s,sp,k,kp,adjr2] = corr_check(x,y)
            x = x(:); y = y(:);
            [p,pp]=corr(x,y,'type','Pearson');
            [s,sp]=corr(x,y,'type','Spearman');
            [k,kp]=corr(x,y,'type','Kendall');
            [~,adjr2]=Utils.regress_line(x,y);
        end
        
        function b = default_low_bounds
            error('deprecated');
            b = [200 40];
        end
        
    end


end