%%TODO: this class will allow visualization of all sessions,
%by mouse, for whichever quantities desired, of the raw curve
%rather than the parameter fit

classdef MultiSessionVisualizer
    methods(Static)
        function Confusion
            load('confusion_agg_190710-161742_0.mat');
            
            tot_num_trials = sum([res.num_trials]);
            tot_cdiff = sum(cat(3, res.sum_CDiff),3)/20; %dividing by number of reps
            normed_cdiff = tot_cdiff./tot_num_trials;
            figure;
            imagesc(normed_cdiff);
            colorbar;
            colormap(bluewhitered);
            xlabel 'Predicted bin'
            ylabel 'Correct bin'
            axis equal;
            xlim([0.5 40.5]); ylim([0.5 40.5]);
            set(gca, 'XTickLabel', {'10R', '20R', '10L', '20L'});
            set(gca, 'YTickLabel', {'10R', '20R', '10L', '20L'});
            line([20.5 20.5], ylim, 'Color', 'k');
            line(xlim, [20.5 20.5], 'Color', 'k');
            figure_format('boxsize', [0.75 0.85]); box on;
            p_ = get(gcf, 'Position');
            set(gcf, 'Position', [p_(1:2), p_(3)*1.5, p_(4)]);
            Utils.create_svg(gcf, 'figure1_svg', 'confusion_diff_both_dirs');
        end
        
        function [series_fits, mouse_names] = Decoding(make_plots)
            if ~exist('make_plots', 'var')
                make_plots = true;
            end
            
            dbfile = 'decoding_all_sess.db';
            conn = sqlite(dbfile);
            samp_size = 20;
            bc = @DecodeTensor.build_command_sess;
            [sess, mouse_names] = DecodeTensor.filt_sess_id_list;
            q = @Utils.cf_p;
            res = q(1,@(s)conn.fetch(bc(s, 'unshuffled', 'MSE', [], 'max')), sess);
            n_sizes = q(1,@(r)double(cell2mat(r(:,1))), res);
            imse = q(1,@(r)1./cell2mat(r(:,3)), res);
            [n_sizes, imse] = Utils.cf_p2(1,...
                @(n,i)MultiSessionVisualizer.regroup(n, i, samp_size),...
                n_sizes, imse);
            
            res_s = q(1,@(s)conn.fetch(bc(s, 'shuffled', 'MSE', [], 'max')), sess);
            n_sizes_s = q(1,@(r)double(cell2mat(r(:,1))), res_s);
            imse_s = q(1,@(r)1./cell2mat(r(:,3)), res_s);
            [n_sizes_s, imse_s] = Utils.cf_p2(1,...
                @(n,i)MultiSessionVisualizer.regroup(n, i, samp_size),...
                n_sizes_s, imse_s);
            assert(isequal(n_sizes, n_sizes_s), 'mismatch between unshuffled and shuffled sampling');
            
            series_fits = q(1,@(s)q(2,@(n,m)createFit_infoSaturation(n(:),mean(m)'), n_sizes, s), {imse, imse_s});
            [I0_fit, I0_conf] = Utils.fit_get(series_fits{1}, 'I_0');
            [I0_fit_s, I0_conf_s] = Utils.fit_get(series_fits{2}, 'I_0');
            
            [N_fit, N_conf] = Utils.fit_get(series_fits{1}, 'N');
            [N_fit_s, N_conf_s] = Utils.fit_get(series_fits{2}, 'N');
            
            save decoding_curves_fits.mat sess mouse_names n_sizes imse imse_s series_fits I0_fit I0_conf I0_fit_s I0_conf_s N_fit N_conf N_fit_s N_conf_s
            
            if make_plots %cancelling unnecessary plots
            MultiSessionVisualizer.plot_series(n_sizes, {imse_s, imse}, {'r', 'b'}, mouse_names, 0.18);
            xlabel 'Number of cells'
            ylabel '1/MSE (cm^{-2})'
            multi_figure_format;
            Utils.create_svg(gcf, 'supplements_svg', 'multi_decoding_IMSE_curves');
            
            figure;
            Utils.bns_groupings(I0_fit, I0_fit_s, I0_conf, I0_conf_s, mouse_names, false);
            xlabel 'Mouse index';
            ylabel(sprintf('I_0 fit value\n(cm^{-2}neuron^{-1})'));
            Utils.specific_format('MBNS');
            Utils.fix_exponent(gca, 'Y', 0);
            Utils.create_svg(gcf, 'supplements_svg', 'multi_I0_fit');
            
            figure;
            Utils.bns_groupings(I0_fit, I0_fit_s, I0_conf, I0_conf_s, mouse_names, true);
            ylim([-Inf Inf]);
            ylabel(sprintf('I_0 fit value\n(cm^{-2}neuron^{-1})'));
            figure_format;
            Utils.fix_exponent(gca, 'Y', 0);
            Utils.create_svg(gcf, 'figure1_svg', 'grouped_I0_fit');
            
            figure;
            Utils.bns_groupings(N_fit, N_fit_s, N_conf, N_conf_s, mouse_names, false);
            set(gca, 'YScale', 'log');
            xlabel 'Mouse index';
            ylabel(sprintf('N fit value\n(neuron)'));
            Utils.specific_format('MBNS');
            Utils.create_svg(gcf, 'supplements_svg', 'multi_N_fit');
            
            figure;
            Utils.bns_groupings(N_fit, N_fit_s, N_conf, N_conf_s, mouse_names, true);
            set(gca, 'YScale', 'log');
            %ylim([-Inf Inf]);
            ylabel(sprintf('N fit value\n(neuron)'));
            figure_format;
            Utils.create_svg(gcf, 'figure1_svg', 'grouped_N_fit');
            
            figure;
            [~,~,select_indices] = DecodeTensor.special_sess_id_list;
            MultiSessionVisualizer.plot_single_filtered(n_sizes, {imse_s, imse}, {'r', 'b'}, select_indices);
            xlabel 'Number of cells'
            ylabel '1/MSE (cm^{-2})'
            %figure_format([0.8125 1.5]);
            figure_format([2 2.5]);
            Utils.create_svg(gcf, 'figure1_svg', 'decoding_IMSE_curves_selected');
            
            figure;
            [~,m_,select_indices] = DecodeTensor.special_sess_id_list;
            %%TODO change colorcell
            %l_ = lines; 
            %l_ = cell2mat(Utils.colorscheme);
            %colorcell = mat2cell(l_(1:numel(select_indices),:), ones(1,numel(select_indices)), 3);
            %colorcell = [Utils.cf_(@(x)x, colorcell), colorcell]';
            MultiSessionVisualizer.plot_single_filtered_sesscolor(n_sizes, {imse_s, imse},...
                [DecodeTensor.mcolor(m_)'; DecodeTensor.mcolor(m_)'],...
                select_indices);
            xlabel 'Number of cells'
            ylabel '1/MSE (cm^{-2})'
            %figure_format([0.8125 1.5]);
            figure_format([2 2.5]);
            Utils.create_svg(gcf, 'figure1_svg', 'decoding_IMSE_curves_selected_colored');
            end %if true
            
            if false
            figure;
            [mouse_identity, agg_n_values, agg_imse_values] =...
                MultiSessionVisualizer.aggregate_sess_per_mouse(n_sizes, imse, mouse_names);
            [mouse_identity_s, agg_n_values_s, agg_imse_values_s] =...
                MultiSessionVisualizer.aggregate_sess_per_mouse(n_sizes, imse_s, mouse_names);
            MultiSessionVisualizer.plot_single_agg(agg_n_values, {agg_imse_values_s, agg_imse_values}, {'r', 'b'}, mouse_identity);
            %almenaux
            xlabel 'Number of cells'
            ylabel '1/MSE (cm^{-2})'
            figure_format([0.8125 1.5].*[2 1.5]);
            Utils.create_svg(gcf, 'figure1_svg', 'decoding_IMSE_curves_mouse_aggregated');
            keyboard;
            agg_imse_mean = Utils.cf_(@(x)cellfun(@mean,x),agg_imse_values);
            agg_imse_mean_s = Utils.cf_(@(x)cellfun(@mean,x),agg_imse_values_s);
            series_fits_agg = q(1,@(s)q(2,@(n,m)createFit_infoSaturation(n(:),m'), agg_n_values, s),...
                {agg_imse_mean, agg_imse_mean_s});
            [agg_I0_fit, agg_I0_conf] = Utils.fit_get(series_fits_agg{1}, 'I_0');
            [agg_I0_fit_s, agg_I0_conf_s] = Utils.fit_get(series_fits_agg{2}, 'I_0');
            
            [agg_N_fit, agg_N_conf] = Utils.fit_get(series_fits_agg{1}, 'N');
            [agg_N_fit_s, agg_N_conf_s] = Utils.fit_get(series_fits_agg{2}, 'N');
            
            figure;
            ballnstick('Unshuffled', 'Shuffled', agg_I0_fit, agg_I0_fit_s, agg_I0_conf, agg_I0_conf_s, 'coloring', DecodeTensor.mcolor(mouse_identity)');
            ylim([-Inf Inf]);
            ylabel(sprintf('I_0 fit value\n(cm^{-2}neuron^{-1})'));
            figure_format;
            Utils.fix_exponent(gca, 'Y', 1);
            Utils.create_svg(gcf, 'figure1_svg', 'agg_I0_fit');
            
            figure;
            ballnstick('Unshuffled', 'Shuffled', agg_N_fit, agg_N_fit_s, agg_N_conf, agg_N_conf_s, 'coloring', DecodeTensor.mcolor(mouse_identity)');
            set(gca, 'YScale', 'log');
            ylim([-Inf Inf]);
            set(gca, 'YTick', [1e2 1e4 1e6]);
            ylabel(sprintf('N fit value\n(neuron)'));
            figure_format;
            Utils.create_svg(gcf, 'figure1_svg', 'agg_N_fit');
            end
        end
        
        function [n, err] = regroup(n_samp, err_samp, samp_size)
            n = unique(n_samp);
            err = zeros(samp_size, numel(n));
            for i = 1:numel(n)
                my_n = n(i);
                my_es = err_samp(n_samp == my_n);
                err(:,i) = my_es(randperm(numel(my_es)) <= samp_size);
            end
        end
        
        function MedLoad
            load('MedLoad_agg_190705-171806_0.mat');
            n_sizes = {res.n_sizes};
            series = {{res.median_loadings}, {res.median_loadings_s}};
            
            series = Utils.cf_(@(m)Utils.cf_(@(x)max(x,[],3),m), series);
            mouse_name = {res.mouse_name};
            MultiSessionVisualizer.plot_series(n_sizes, series, {'b','r'}, mouse_name);
            axs = findall(gcf, 'type', 'axes');
            set(axs, 'YScale', 'log');
            set(axs, 'XScale', 'log');
            xlabel 'Number of cells'
            ylabel 'max_i|cos(PC_i, \Delta\mu)|'
            multi_figure_format;
            Utils.create_svg(gcf, 'supplements_svg', 'multi_medload');
            
            n_c = 50;
            fit_func = @(x,y)fit(log10(x(x>=n_c))',log10(mean(y(:,x>=n_c)))', 'poly1');
            fr_ = Utils.cf_(fit_func, n_sizes, series{1});
            fr_s = Utils.cf_(fit_func, n_sizes, series{2});
            figure;
            [rate_f, rate_f_conf] = Utils.fit_get(fr_, 'p1');
            [rate_f_s, rate_f_s_conf] = Utils.fit_get(fr_s, 'p1');
            Utils.bns_groupings(rate_f, rate_f_s, rate_f_conf, rate_f_s_conf, mouse_name, false);
            xlabel 'Mouse index'
            ylabel 'Fit exponent'
            Utils.specific_format('MBNS');
            Utils.create_svg(gcf, 'supplements_svg', 'multi_medload_exponents');
            
            figure;
            Utils.bns_groupings(rate_f, rate_f_s, rate_f_conf, rate_f_s_conf, mouse_name, true);
            hold on;
            line(xlim-0.5, [0 0], 'Color', 'k', 'LineStyle', '-');
            line(xlim, [-0.5 -0.5], 'Color', 'k', 'LineStyle', ':');
            ylabel 'Fit exponent'
            ylim([-Inf Inf]);
            set(gca, 'XTickLabels', {'Unsh.', 'Sh.'});
            %figure_format([0.8 1]/2, 'fontsize', 4);
            Utils.specific_format('inset');
            Utils.create_svg(gcf, 'figure2_svg', 'group_medload_exponents');
            
            figure;
            %boxplot([rate_f(:), rate_f_s(:)], {'Unsh.', 'Sh.'});
            Utils.basic_boxplot('Unsh.','Sh.',rate_f,rate_f_s);
            hold on;
            line(xlim-0.5, [0 0], 'Color', 'k', 'LineStyle', '-');
            line(xlim, [-0.5 -0.5], 'Color', 'k', 'LineStyle', ':');
            ylim([-Inf Inf]);
            %ylabel('Slope past 50 cells');
            figure_format([0.8 1]/2, 'fontsize', 5);
            Utils.create_svg(gcf, 'figure2_svg', 'medload_selected_inset');
            
            figure;
            [~,m_,select_indices] = DecodeTensor.special_sess_id_list;
            MultiSessionVisualizer.plot_single_filtered(n_sizes, series, {'b','r'}, select_indices);
            set(gca, 'YScale', 'log');
            set(gca, 'XScale', 'log');
            xlabel 'Number of cells'
            ylabel 'max_i|cos(PC_i, Dm)|'
            xlim([1 500]);
            ylim([-Inf 1]);
            figure_format([1 1.4]);
            %figure_format;
            Utils.create_svg(gcf, 'figure2_svg', 'medload_selected');
            
            figure;
            MultiSessionVisualizer.plot_single_filtered_sesscolor(n_sizes,...
                series, [DecodeTensor.mcolor(m_)'; DecodeTensor.mcolor(m_)'], select_indices);
            set(gca, 'YScale', 'log');
            set(gca, 'XScale', 'log');
            xlabel 'Number of cells'
            ylabel 'max_i|cos(PC_i, Dm)|'
            xlim([1 500]);
            ylim([-Inf 1]);
            figure_format([1 1.4]);
            %figure_format;
            Utils.create_svg(gcf, 'figure2_svg', 'medload_selected_colored');
        end
        
        function SnN(normed, filt_num)
            %load('signal_and_noise_final.mat');
            load('SnN1000_agg_190708-235729_0.mat');
            [~, mouse_list] = DecodeTensor.filt_sess_id_list;
            
            dm2_full = {res.dm2}; sm_full = {res.sm}; sms_full = {res.sms}; n_sizes_full = {res.n_sizes};
            series = {dm2_full, sm_full, sms_full};
            series_colors = {'k', 'b', 'r'};
            
            if normed
                vs = Utils.cf_(@(n,m)Utils.fitaline(n,m), n_sizes_full, dm2_full);
                %sm_full = cellfun(@(s,v)s./v, sm_full, vs, 'UniformOutput', false);
                %sms_full = cellfun(@(s,v)s./v, sms_full, vs, 'UniformOutput', false);
                %dm2_full = cellfun(@(s,v)s./v, dm2_full, vs, 'UniformOutput', false);
                func = @(x) Utils.cf_(@(s,v)s./v, x, vs);
                series = Utils.cf_(func, series);
            end
            
            if filt_num
                f_ = cellfun(@(n)n(end)>=200, n_sizes_full);
                n_sizes_full = n_sizes_full(f_);
                %dm2_full = dm2_full(f_);
                %sms_full = sms_full(f_);
                %sm_full = sm_full(f_);
                series = Utils.cf_(@(x)x(f_), series);
                mouse_list = mouse_list(f_);
            end
            
            MultiSessionVisualizer.plot_series(n_sizes_full,...
                series, series_colors, mouse_list);
            xlabel 'Number of cells'
            if normed
                ylabel(sprintf('Distance^2\n(in units of cells)'));
            else
                ylabel(sprintf('Distance^2\n(on \\DeltaF/F values)'));
            end
            multi_figure_format;
            Utils.create_svg(gcf, 'supplements_svg', 'multi_signal_and_noise_growth');
            
            q=@Utils.cf_p;
            n_c = 100;
            series_fits = q(1,@(s)q(2,@(n,m)fit(n(n>=n_c)',mean(m(:,n>=n_c))','poly1'), n_sizes_full, s), series);
            %progressbar('series', 'fittings');
            %series_fits = q(1, @(s)q(2, @(n,m)Utils.fit_slopechanger(n, mean(m)), n_sizes_full, s), series);
            [sm_slope, sm_conf] = Utils.fit_get(series_fits{2}, 'p1');%'r_f');
            [sms_slope, sms_conf] = Utils.fit_get(series_fits{3}, 'p1');%'r_f');
            [sms_intercept, sms_intercept_conf] = Utils.fit_get(series_fits{3}, 'p2');
            %[sm_slope, sm_conf] = Utils.fit_get(series_fits{2}, 'r_f');%'p1');
            %[sms_slope, sms_conf] = Utils.fit_get(series_fits{3}, 'r_f');%'p1');
            figure;
            Utils.bns_groupings(sm_slope, sms_slope, sm_conf, sms_conf, mouse_list, false);
            xlabel 'Mouse index';
            ylabel(sprintf('s^2 along Dm\nrate of change'));
            Utils.specific_format('MBNS');
            Utils.create_svg(gcf, 'supplements_svg', 'multi_noise_rate_of_change');
            
            figure;
            Utils.bns_groupings(sm_slope, sms_slope, sm_conf, sms_conf, mouse_list, true);
            ylim([-Inf Inf]);
            ylabel(sprintf('s^2 along Dm\nrate of change'));
            figure_format('factor', 1.6);
            Utils.create_svg(gcf, 'figure2_svg', 'grouped_noise_rate_of_change');
            
            figure;
            Utils.basic_doublehist('Unshuffled', 'Shuffled', sm_slope, sms_slope, -0.1:0.05:1.1);
            xlim([-0.1 1.1]);
            xlabel(sprintf('s^2 along Dm\nrate of change'));
            figure_format('factor', 1.3);
            Utils.create_svg(gcf, 'figure2_svg', 'hist_noise_rate_of_change');
            
            %pause;
            
            figure;
            Utils.horiz_boxplot('Unsh.', 'Sh.', sm_slope, sms_slope);
            %l_ = refline(0, 1); l_.Color = 'g'; l_.LineStyle = '-';
            %l_ = refline(0, 0); l_.Color = 'g'; l_.LineStyle = '-';
            %ylabel(sprintf('\\sigma^2 along \\Delta\\mu\nrate of change'));
            line([-0.1 -0.1], ylim, 'Color', 'g');
            line([1.1 1.1], ylim, 'Color', 'g');
            Utils.specific_format('inset');
            Utils.create_svg(gcf, 'figure2_svg', 'boxplot_noise_rate_of_change');
            
            figure;
            [~,my_mice,select_indices] = DecodeTensor.special_sess_id_list;
            if filt_num
                [~, my_mice] = DecodeTensor.filt_sess_id_list;
                select_indices_filt = false(1,numel(my_mice));
                select_indices_filt(select_indices) = true;
                select_indices_filt = select_indices_filt(f_);
                select_indices = find(select_indices_filt);
                my_mice = my_mice(select_indices_filt);
            end
            MultiSessionVisualizer.plot_single_filtered(n_sizes_full, series([1 3 2]), series_colors([1 3 2]), select_indices);
            xlabel 'Number of cells'
            if normed
                ylabel(sprintf('Distance^2\n(in units of cells)'));
            else
                ylabel(sprintf('Distance^2\n(on \\DeltaF/F values)'));
            end
            text(20, 370, '(Dm)^2', 'Color', 'k', 'HorizontalAlignment', 'left');
            text(20, 470, 's^2 along Dm', 'Color', 'b', 'HorizontalAlignment', 'left');
            text(20, 570, 's^2 along Dm (Shuffled)', 'Color', 'r', 'HorizontalAlignment', 'left');
            %figure_format('factor', 1.3);
            figure_format([1 1.4], 'factor', 1.2);
            Utils.create_svg(gcf, 'figure2_svg', 'signal_and_noise_growth');
            
            figure;
            assert(numel(select_indices) == numel(DecodeTensor.mcolor(my_mice)), 'should have 12 mice');
            MultiSessionVisualizer.plot_single_filtered_sesscolor(n_sizes_full,...
                series([1 3 2]), [repmat({'k'}, 1, numel(select_indices));...
                DecodeTensor.mcolor(my_mice)';...
                DecodeTensor.mcolor(my_mice)'], select_indices);
            xlabel 'Number of cells'
            if normed
                ylabel(sprintf('Distance^2\n(in units of cells)'));
            else
                ylabel(sprintf('Distance^2\n(on \\DeltaF/F values)'));
            end
            text(20, 370, '(Dm)^2', 'Color', 'k', 'HorizontalAlignment', 'left');
            text(20, 470, 's^2 along Dm', 'Color', 'b', 'HorizontalAlignment', 'left');
            text(20, 570, 's^2 along Dm (Shuffled)', 'Color', 'r', 'HorizontalAlignment', 'left');
            %figure_format('factor', 1.3);
            figure_format([1 1.4], 'factor', 1.2);
            Utils.create_svg(gcf, 'figure2_svg', 'signal_and_noise_growth_colored');
            %return;%TODO remove
            %load('fit_result_record.mat');
            %bagopo
            dotsize = 4;
            [series_fits, mouse_names] = MultiSessionVisualizer.Decoding(false);
            
            fitresult = series_fits{1};
            fitresult_s = series_fits{2};
            if filt_num
                fitresult = fitresult(f_);
                fitresult_s = fitresult_s(f_);
                mouse_names = mouse_names(f_);
                %fitresult_d = fitresult_d(f_);
            end
            [I0_fit_value, I0_upper] = Utils.fit_get(fitresult, 'I_0');
            [I0_fit_value_s, I0_upper_s] = Utils.fit_get(fitresult_s, 'I_0');
            [N_fit_value, N_upper] = Utils.fit_get(fitresult, 'N');
            
            good_fit_filter = (I0_upper < I0_fit_value) &...
                (I0_upper_s < I0_fit_value_s) &...
                (N_upper < N_fit_value);
            g_ = good_fit_filter;
            disp(find(~g_));
            figure;
            %scatter(I0_fit_value.*N_fit_value, 1./sm_slope, 4, 'b');
            limit_uncertainty = sqrt((I0_fit_value.*(N_upper)).^2 + (N_fit_value.*(I0_upper)).^2);
            inv_sm_slope_uncertainty = sm_conf ./ sm_slope.^2;
            hold on;
            
            errorbar(I0_fit_value(g_).*N_fit_value(g_), 1./sm_slope(g_), inv_sm_slope_uncertainty(g_), inv_sm_slope_uncertainty(g_),...
                limit_uncertainty(g_), limit_uncertainty(g_), 'LineStyle', 'none', 'Color', 'k', 'CapSize', 1);
            scatter(I0_fit_value(g_).*N_fit_value(g_), 1./sm_slope(g_), dotsize, DecodeTensor.mcolor(mouse_names(g_), false), 'filled');
            [fitresult, adjr2] = Utils.regress_line(I0_fit_value(g_).*N_fit_value(g_), 1./sm_slope(g_));
            h_ = plot(fitresult); legend off
            h_.Color = 'b';
            xlim([-Inf 0.15]);
            text(0.1, 1, sprintf('adj. R^2 = %.2f', adjr2));
            xlabel 'IMSE limit I_0N';
            ylabel(sprintf('Inverse s^2\nrate of change'));
            fprintf('imse limit regression, N = %d\n', numel(I0_fit_value(g_)));
            figure_format('factor', 1.6);
            Utils.create_svg(gcf, 'figure2_svg', 'imse_limit_regression');
            
            
            figure;
            %scatter(I0_fit_value_s, 1./sms_intercept, 'r');
            hold on;
            inv_intercept_errb = sms_intercept_conf./sms_intercept.^2;
            
            errorbar(I0_fit_value_s(g_), 1./sms_intercept(g_), inv_intercept_errb(g_), inv_intercept_errb(g_), I0_upper_s(g_), I0_upper_s(g_), 'LineStyle', 'none', 'Color', 'k', 'CapSize', 1);
            scatter(I0_fit_value_s(g_), 1./sms_intercept(g_), dotsize, DecodeTensor.mcolor(mouse_names(g_), false), 'filled');
            [fitresult, adjr2] = Utils.regress_line(I0_fit_value_s(g_), 1./sms_intercept(g_));
            plot(fitresult); legend off
            text(7e-4, 0.055, sprintf('adj. R^2 = %.2f', adjr2));
            xlabel 'I_0 fit value'
            ylabel 'Asymptotic 1/s^2'
            set(gca, 'XTickLabel', arrayfun(@Utils.my_fmt, get(gca, 'XTick') ,'UniformOutput', false));
            fprintf('I0 value regression, N = %d\n', numel(I0_fit_value_s(g_)));
            figure_format('factor', 1.6);
            Utils.create_svg(gcf, 'figure2_svg', 'I0_value_regression');
        end
        
        function plot_series(n_sizes, series_cell, color_cell, mouse_list, max_y_val)
            if ~exist('max_y_val', 'var')
                max_y_val = Inf;
            end
            if numel(max_y_val) > 1
                max_y_val = max_y_val(end);
            end
            figure;
            mouse_names = unique(mouse_list);
            num_mice = numel(mouse_names);
            n_rows = round(sqrt(num_mice));
            n_cols = ceil(num_mice/n_rows);
            for m_i = 1:num_mice
                subplot(n_rows, n_cols, m_i);
                f_ = strcmp(mouse_names{m_i}, mouse_list);
                n_ = n_sizes(f_);
                for j = 1:numel(series_cell)
                    s_ = series_cell{j}(f_);
                    c_ = color_cell{j};
                    for k = 1:numel(s_)
                        Utils.neuseries(n_{k}, s_{k}, c_);
                        hold on;
                    end %session in mouse
                end %quantity shown
                xlim([0 500]);
                ylim([0 max_y_val]);
                title(mouse_names{m_i});
            end %mice id
        end %func
        
        function plot_single_filtered(n_sizes, series_cell, color_cell, filter_selection)
            n_sizes = n_sizes(filter_selection);
            series_cell = Utils.cf_(@(x)x(filter_selection), series_cell);
            for j = 1:numel(series_cell)
                s_ = series_cell{j};
                c_ = color_cell{j};
                for k = 1:numel(s_)
                    Utils.neuseries(n_sizes{k}, s_{k}, c_);
                    hold on;
                end %session
            end %quantity shown
            xlim([0 500]);
            ylim([0 Inf]);
        end%func
        
        function plot_single_agg(n_sizes, series_cell, color_cell, mouse_identity)
            for j = 1:numel(series_cell)
                s_ = series_cell{j};
                c_ = color_cell{j};
                for k = 1:numel(s_)
                    s_set = s_{k};
                    s_mean = cellfun(@mean, s_set);
                    s_err = 1.96.*cellfun(@(x)std(x)./sqrt(length(x)), s_set);
                    hold on;
                    shadedErrorBar(n_sizes{k}, s_mean, s_err, 'lineprops', c_);
                end
            end
        end
        
        function plot_single_filtered_sesscolor(n_sizes, series_cell, color_cell, filter_selection)
            n_sizes = n_sizes(filter_selection);
            series_cell = Utils.cf_(@(x)x(filter_selection), series_cell);
            for j = 1:numel(series_cell)
                s_ = series_cell{j};
                %c_ = color_cell{j};
                for k = 1:numel(s_)
                    h_ = Utils.neuseries(n_sizes{k}, s_{k}, 'k');
                    set(h_.edge, 'Color', color_cell{j,k});
                    h_.patch.FaceColor = color_cell{j,k};
                    h_.patch.EdgeColor = color_cell{j,k};
                    h_.mainLine.Color = color_cell{j,k};
                    hold on;
                end %session
            end %quantity shown
            xlim([0 500]);
            ylim([0 Inf]);
        end%func
        
        function [my_mice, n_vals, corresponding_mean_ms] = aggregate_sess_per_mouse(ns, ms, mouse_name)
            ns = ns(:); ms = ms(:); mouse_name = mouse_name(:);
            my_mice = unique(mouse_name);
            for m_i = 1:numel(my_mice)
                my_mouse = my_mice{m_i};
                f_ = strcmp(mouse_name, my_mouse);
                my_ns = ns(f_);
                my_ms = ms(f_);
                n_vals{m_i} = unique(cell2mat(Utils.cf_(@(x)x(:)', my_ns(:)')));
                for n_i = 1:numel(n_vals{m_i})
                    my_n = n_vals{m_i}(n_i);
                    corresponding_mean_ms{m_i}{n_i} = ...
                        cell2mat(Utils.cf_(@(f, m_) mean(m_(:,f)),...
                        Utils.cf_(@(n_) n_ == my_n, my_ns), my_ms)');
                end
            end
        end
    end%methods
end%classdef
