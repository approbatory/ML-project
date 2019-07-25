classdef PanelGenerator
    methods(Static)
        function confusion(remake)
            if ~exist('remake', 'var')
                remake = false;
            end
            savedir = 'figure1_pdf/confusion';
            if ~exist(savedir, 'dir')
                mkdir(savedir);
            end
            
            ap = @(x) fullfile(savedir, x);
            
            load('confusion_single_session_agg_190723-125628_0.mat');
            nt = res(1).num_trials;
            C_pct = cat(3, res.C)/nt*100;
            C_s_pct = cat(3, res.C_s)/nt*100;
            nb = size(C_pct, 1);
            
            fname = ap('confusion_unshuffled.pdf');
            if remake || ~exist(fname, 'file')
                figure('FileName', fname);
                imagesc(squeeze(mean(C_pct,3)), [0 100]);
                xlabel 'Predicted bin'
                ylabel 'Correct bin'
                axis equal;
                xlim([1 nb] + [-0.5 0.5]);
                set(gca, 'XTickLabel', []);%{'60 cm', '120 cm', '60 cm', '120 cm'});
                set(gca, 'YTickLabel', []);%{'60 cm', '120 cm', '60 cm', '120 cm'});
                line([nb nb]/2+0.5, ylim, 'Color', 'w');
                line(xlim, [nb nb]/2+0.5, 'Color', 'w');
                h_ = colorbar; ylabel(h_, 'Confusion (%)', 'Rotation', 270);
                Utils.specific_format('confusion');
                Utils.printto;
            end
            
            fname = ap('confusion_shuffled.pdf');
            if remake || ~exist(fname, 'file')
                figure('FileName', fname);
                imagesc(squeeze(mean(C_s_pct,3)), [0 100]);
                xlabel 'Predicted bin'
                ylabel 'Correct bin'
                axis equal;
                xlim([1 nb] + [-0.5 0.5]);
                set(gca, 'XTickLabel', []);%{'60 cm', '120 cm', '60 cm', '120 cm'});
                set(gca, 'YTickLabel', []);%{'60 cm', '120 cm', '60 cm', '120 cm'});
                line([nb nb]/2+0.5, ylim, 'Color', 'w');
                line(xlim, [nb nb]/2+0.5, 'Color', 'w');
                h_ = colorbar; ylabel(h_, 'Confusion (%)', 'Rotation', 270);
                Utils.specific_format('confusion');
                Utils.printto;
            end
            
            fname = ap('confusion_diff.pdf');
            if remake || ~exist(fname, 'file')
                figure('FileName', fname);
                imagesc(squeeze(mean(C_s_pct,3)) - squeeze(mean(C_pct,3)));
                xlabel 'Predicted bin'
                ylabel 'Correct bin'
                axis equal;
                xlim([1 nb] + [-0.5 0.5]);
                set(gca, 'XTickLabel', []);%{'60 cm', '120 cm', '60 cm', '120 cm'});
                set(gca, 'YTickLabel', []);%{'60 cm', '120 cm', '60 cm', '120 cm'});
                line([nb nb]/2+0.5, ylim, 'Color', 'k');
                line(xlim, [nb nb]/2+0.5, 'Color', 'k');
                colormap(bluewhitered);
                h_ = colorbar; ylabel(h_, 'Confusion difference (%)', 'Rotation', 270);
                Utils.specific_format('confusion');
                Utils.printto;
            end
            
            fname = ap('bin_error_rate.pdf');
            if remake || ~exist(fname, 'file') %TODO this should have errorbars
                figure('FileName', fname);
                m = 100-diag(squeeze(mean(C_pct,3)));
                e = diag(squeeze(std(C_pct,[],3)))./sqrt(size(C_pct,3));
                errorbar(m, e, 'b', 'CapSize', 1);
                xlabel 'Place bin'
                ylabel 'Error rate (%)'
                hold on
                m_s = 100-diag(squeeze(mean(C_s_pct,3)));
                e_s = diag(squeeze(std(C_s_pct,[],3)))./sqrt(size(C_s_pct,3));
                errorbar(m_s, e_s, 'r', 'CapSize', 1);
                errorbar(m - m_s, sqrt(e.^2 + e_s.^2), 'k', 'DisplayName', 'Difference', 'CapSize', 1);
                line([20 20], ylim, 'Color', 'k');
                figure_format([0.85 0.85]);
                Utils.printto;
            end
        end
        
        function decoding_curves(varargin)
            p = inputParser;
            p.addOptional('remake', false, @islogical);
            p.addOptional('recompute', false, @islogical);
            
            p.parse(varargin{:});
            
            savedir = 'figure1_pdf/decoding_curves';
            if ~exist(savedir, 'dir')
                mkdir(savedir);
            end
            
            ap = @(x) fullfile(savedir, x);
            
            
            save_file = 'decoding_curves_fits.mat';
            if p.Results.recompute || ~exist(save_file, 'file')
                dbfile = 'decoding_all_sess.db';
                conn = sqlite(dbfile);
                samp_size = 80;
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
                
                %series_fits = q(1,@(s)q(2,@(n,m)createFit_infoSaturation(n(:),mean(m)'), n_sizes, s), {imse, imse_s});
                
                [series_fits{1}, series_gof{1}] = Utils.cf_p2(2,@(n,m)createFit_infoSaturation(n(:),mean(m)'), n_sizes, imse);
                progressbar(0.5/2);
                [series_fits{2}, series_gof{2}] = Utils.cf_p2(2,@(n,m)createFit_infoSaturation(n(:),mean(m)'), n_sizes, imse_s);
                progressbar(1/2);
                [series_exp_fits{1}, series_exp_gof{1}] = Utils.cf_p2(2,@(n,m)createFit_exp(n(:),mean(m)'), n_sizes, imse);
                progressbar(1.5/2);
                [series_exp_fits{2}, series_exp_gof{2}] = Utils.cf_p2(2,@(n,m)createFit_exp(n(:),mean(m)'), n_sizes, imse_s);
                progressbar(2/2);
                [I0_fit, I0_conf] = Utils.fit_get(series_fits{1}, 'I_0');
                [I0_fit_s, I0_conf_s] = Utils.fit_get(series_fits{2}, 'I_0');
                
                [N_fit, N_conf] = Utils.fit_get(series_fits{1}, 'N');
                [N_fit_s, N_conf_s] = Utils.fit_get(series_fits{2}, 'N');
                
                save(save_file, 'sess', 'mouse_names', 'n_sizes',...
                    'imse', 'imse_s', 'series_fits', 'series_gof', 'I0_fit', 'I0_conf',...
                    'I0_fit_s', 'I0_conf_s', 'N_fit', 'N_conf',...
                    'N_fit_s', 'N_conf_s', 'series_exp_fits', 'series_exp_gof');
            else
                load(save_file);
            end
            
            r2 = cellfun(@(x)x.rsquare, series_gof{1});
            r2_s = cellfun(@(x)x.rsquare, series_gof{2});
            
            r2_exp = cellfun(@(x)x.rsquare, series_exp_gof{1});
            r2_s_exp = cellfun(@(x)x.rsquare, series_exp_gof{2});
            
            fprintf('Unshuf R^2: %f - %f, median: %f\n', min(r2), max(r2), median(r2));
            fprintf('Shuf R^2: %f - %f, median: %f\n', min(r2_s), max(r2_s), median(r2_s));
            fprintf('Unshuf (exp) R^2: %f - %f, median: %f\n', min(r2_exp), max(r2_exp), median(r2_exp));
            fprintf('Shuf (exp) R^2: %f - %f, median: %f\n', min(r2_s_exp), max(r2_s_exp), median(r2_s_exp));
            
            
            fname = ap('decoding_curve_fit.pdf');
            if p.Results.remake || ~exist(fname, 'file')
                figure('FileName', fname);
                show_mice = {'Mouse2022', 'Mouse2024', 'Mouse2028'};
                
                [~,m_,sp_] = DecodeTensor.special_sess_id_list;
                show_filter = ismember(m_, show_mice);
                sp_ = sp_(show_filter);
                %MultiSessionVisualizer.plot_single_filtered(n_sizes, {imse_s, imse}, {'r', 'b'}, find(show_filter));
                for j = 1:numel(sess)
                    if ismember(j, sp_)
                        n = n_sizes{j};
                        %i = imse{i};
                        i_s = imse_s{j};
                        %plot(series_fits{2}{j}, 'r');
                        n_f = 1:500;
                        plot(n_f, I0_fit_s(j).*n_f./(1 + n_f./N_fit_s(j)), 'r');
                        hold on;
                        errorbar(n, mean(i_s), std(i_s)./sqrt(size(i_s,1)),...
                            'r', 'Capsize', 1, 'LineStyle', 'none');
                    end
                end
                for j = 1:numel(sess)
                    if ismember(j, sp_)
                        n = n_sizes{j};
                        i = imse{j};
                        %plot(series_fits{1}{j}, 'b');
                        n_f = 1:500;
                        plot(n_f, I0_fit(j).*n_f./(1 + n_f./N_fit(j)), 'b');
                        hold on;
                        errorbar(n, mean(i), std(i)./sqrt(size(i,1)), 'b',...
                            'Capsize', 1, 'LineStyle', 'none');
                    end
                end
                xlabel 'Number of cells'
                ylabel '1/MSE (cm^{-2})'
                ylim([0 0.16]);
                legend off
                %figure_format([0.8125 1.5]);
                figure_format([2 2.5]);
                Utils.printto;
                
                [pref, fn, ext] = fileparts(fname);
                figure('FileName', fullfile(pref, ['inset' fn ext]));
                for j = 1:numel(sess)
                    if ismember(j, sp_)
                        n = n_sizes{j};
                        %i = imse{i};
                        i_s = imse_s{j};
                        %plot(series_fits{2}{j}, 'r');
                        n_f = 1:500;
                        plot(n_f, (I0_fit_s(j).*n_f./(1 + n_f./N_fit_s(j))).^(-1/2), 'r');
                        hold on;
                        errorbar(n, mean(i_s.^(-1/2)), std(i_s.^(-1/2))./sqrt(size(i_s,1)),...
                            'r', 'Capsize', 0.4, 'LineStyle', 'none');
                    end
                end
                for j = 1:numel(sess)
                    if ismember(j, sp_)
                        n = n_sizes{j};
                        i = imse{j};
                        %plot(series_fits{1}{j}, 'b');
                        n_f = 1:500;
                        plot(n_f, (I0_fit(j).*n_f./(1 + n_f./N_fit(j))).^(-1/2), 'b');
                        hold on;
                        errorbar(n, mean(i.^(-1/2)), std(i.^(-1/2))./sqrt(size(i,1)), 'b',...
                            'Capsize', 0.4, 'LineStyle', 'none');
                    end
                end
                xlabel 'Number of cells'
                ylabel 'RMS Error (cm)'
                ylim([2 50]);
                set(gca, 'YScale', 'log');
                set(gca, 'YTick', [1 2 5 10 20 50]);
                legend off
                figure_format('boxsize', [0.6 0.8], 'fontsize', 5);
                Utils.printto;
            end
            
            fname = ap('grouped_I0_fit.pdf');
            if p.Results.remake || ~exist(fname, 'file')
                figure('FileName', fname);
                Utils.bns_groupings(I0_fit, I0_fit_s, I0_conf, I0_conf_s, mouse_names, true);
                ylim([-Inf Inf]);
                ylabel(sprintf('I_0 fit value\n(cm^{-2}neuron^{-1})'));
                figure_format;
                Utils.fix_exponent(gca, 'Y', 0);
                Utils.printto;
            end
            
            fname = ap('grouped_N_fit.pdf');
            if p.Results.remake || ~exist(fname, 'file')
                figure('FileName', fname);
                Utils.bns_groupings(N_fit, N_fit_s, N_conf, N_conf_s, mouse_names, true);
                set(gca, 'YScale', 'log');
                %ylim([-Inf Inf]);
                ylabel(sprintf('N fit value\n(neuron)'));
                figure_format;
                Utils.printto;
            end
        end
        
        function medload
            load('MedLoad_agg_190705-171806_0.mat');
            n_sizes = {res.n_sizes};
            series = {{res.median_loadings}, {res.median_loadings_s}};
            
            series = Utils.cf_(@(m)Utils.cf_(@(x)max(x,[],3),m), series);
            mouse_name = {res.mouse_name};
            
            show_mice = {'Mouse2022', 'Mouse2024', 'Mouse2028'};
            
            [~,m_,sp_] = DecodeTensor.special_sess_id_list;
            show_filter = ismember(m_, show_mice);
            sp_ = sp_(show_filter);
            figure('FileName', 'figure2_pdf/medload/medload_curves.pdf');
            MultiSessionVisualizer.plot_single_filtered(n_sizes, series, {'b', 'r'}, sp_);
            set(gca, 'XScale', 'log');
            set(gca, 'YScale', 'log');
            xlabel 'Number of cells'
            ylabel 'max_i|cos(PC_i, Dm)|'
            xlim([1 500]);
            ylim([-Inf 1]);
            figure_format([1 1.4]);
            Utils.printto;
            
            n_c = 50;
            fit_func = @(x,y)fit(log10(x(x>=n_c))',log10(mean(y(:,x>=n_c)))', 'poly1');
            fr_ = Utils.cf_(fit_func, n_sizes, series{1});
            fr_s = Utils.cf_(fit_func, n_sizes, series{2});
            
            [rate_f, rate_f_conf] = Utils.fit_get(fr_, 'p1');
            [rate_f_s, rate_f_s_conf] = Utils.fit_get(fr_s, 'p1');
            figure('FileName', 'figure2_pdf/medload/inset.pdf');
            Utils.bns_groupings(rate_f, rate_f_s, rate_f_conf, rate_f_s_conf, mouse_name, true);
            hold on;
            line(xlim-0.5, [0 0], 'Color', 'k', 'LineStyle', '-');
            line(xlim, [-0.5 -0.5], 'Color', 'k', 'LineStyle', ':');
            ylabel 'Fit exponent'
            ylim([-Inf Inf]);
            set(gca, 'XTickLabels', {'Unsh.', 'Sh.'});
            %figure_format([0.8 1]/2, 'fontsize', 4);
            Utils.specific_format('inset');
            Utils.printto;
        end
        
        function signal_and_noise
            normed = true;
            filt_num = false;
            
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
            
             q=@Utils.cf_p;
            n_c = 100;
            series_fits = q(1,@(s)q(2,@(n,m)fit(n(n>=n_c)',mean(m(:,n>=n_c))','poly1'), n_sizes_full, s), series);
            %progressbar('series', 'fittings');
            %series_fits = q(1, @(s)q(2, @(n,m)Utils.fit_slopechanger(n, mean(m)), n_sizes_full, s), series);
            [sm_slope, sm_conf] = Utils.fit_get(series_fits{2}, 'p1');%'r_f');
            [sms_slope, sms_conf] = Utils.fit_get(series_fits{3}, 'p1');%'r_f');
            [sms_intercept, sms_intercept_conf] = Utils.fit_get(series_fits{3}, 'p2');
            
            figure('FileName', 'figure2_pdf/signal_and_noise/noise_rate_of_change.pdf');
            Utils.bns_groupings(sm_slope, sms_slope, sm_conf, sms_conf, mouse_list, true);
            ylim([-Inf Inf]);
            ylabel(sprintf('s^2 along Dm\nrate of change'));
            figure_format('factor', 1.6);
            Utils.printto;
            %Utils.create_svg(gcf, 'figure2_svg', 'grouped_noise_rate_of_change');
            
            figure('FileName', 'figure2_pdf/signal_and_noise/noise_curves.pdf');
            show_mice = {'Mouse2022', 'Mouse2024', 'Mouse2028'};
            
            [~,m_,sp_] = DecodeTensor.special_sess_id_list;
            show_filter = ismember(m_, show_mice);
            sp_ = sp_(show_filter);
            MultiSessionVisualizer.plot_single_filtered(n_sizes_full, series([1 3 2]), series_colors([1 3 2]), sp_);
            xlabel 'Number of cells'
            if normed
                ylabel(sprintf('Noise variance/K\n(in units of cells)'));
            else
                ylabel(sprintf('Noise variance\n(on \\DeltaF/F values)'));
            end
            text(20, 370, '(Dm)^2', 'Color', 'k', 'HorizontalAlignment', 'left');
            text(20, 470, 's^2 along Dm', 'Color', 'b', 'HorizontalAlignment', 'left');
            text(20, 570, 's^2 along Dm (Shuffled)', 'Color', 'r', 'HorizontalAlignment', 'left');
            %figure_format('factor', 1.3);
            figure_format([1 1.4], 'factor', 1.2);
            Utils.printto;
            
            fit_savefile = 'decoding_curves_fits.mat';
            if ~exist(fit_savefile, 'file')
                PanelGenerator.decoding_curves;
            end
            load(fit_savefile);
            dotsize = 4;
            
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
            figure('FileName', 'figure2_pdf/signal_and_noise/imse_limit_regression.pdf');
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
            %Utils.create_svg(gcf, 'figure2_svg', 'imse_limit_regression');
            Utils.printto;
            
            
            figure('FileName', 'figure2_pdf/signal_and_noise/I0_value_regression.pdf');
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
            %Utils.create_svg(gcf, 'figure2_svg', 'I0_value_regression');
            Utils.printto;
        end
    end
end