classdef PanelGenerator
    methods(Static)
        function plot_confusion(C_pct, cbar_label, clim_exists)
            %Plot a confusion matrix, or difference between confusion
            %matrices.
            %Inputs:
            %   C_pct: square confusion matrix where elements are percentages
            %   cbar_label: string to display as a label to the colorbar
            %   clim_exists (optional, default true): whether or not to
            %   show the color range from 0 to 100 rather than min to max
            nb = size(C_pct, 1);
            if ~exist('clim_exists', 'var') || clim_exists
                imagesc(squeeze(mean(C_pct,3)), [0 100]);
            else
                imagesc(squeeze(mean(C_pct,3)));
            end
            xlabel 'Predicted bin'
            ylabel 'Correct bin'
            axis equal;
            xlim([1 nb] + [-0.5 0.5]);
            set(gca, 'XTickLabel', []);%{'60 cm', '120 cm', '60 cm', '120 cm'});
            set(gca, 'YTickLabel', []);%{'60 cm', '120 cm', '60 cm', '120 cm'});
            line([nb nb]/2+0.5, ylim, 'Color', 'w');
            line(xlim, [nb nb]/2+0.5, 'Color', 'w');
            h_ = colorbar; ylabel(h_, cbar_label, 'Rotation', 270);
            Utils.specific_format('confusion');
        end
        
        function [n_sizes, imse] = db_imse_reader(conn, setting, sess, samp_size)
            %Read from the decoding database.
            %Inputs:
            %   conn: a valid connection to a sqlite database
            %   setting: either 'unshuffled', 'shuffled', or 'diagonal'
            %   sess: a string cell of session codes denoting which
            %       sessions to read out
            %   samp_size: how many samples to load, for a given set of
            %       parameters, e.g. 20 or 80
            bc = @DecodeTensor.build_command_sess;
            q = @Utils.cf_p;
            res = q(1,@(s)conn.fetch(bc(s, setting, 'MSE', [], 'max')), sess);
            n_sizes = q(1,@(r)double(cell2mat(r(:,1))), res);
            imse = q(1,@(r)1./cell2mat(r(:,3)), res);
            [n_sizes, imse] = Utils.cf_p2(1,...
                @(n,i)MultiSessionVisualizer.regroup(n, i, samp_size),...
                n_sizes, imse);
        end
        
        function [n_sizes, imse, mask] = db_imse_reader_safe(conn, setting, sess, samp_size)
            %Read from the decoding database.
            %Inputs:
            %   conn: a valid connection to a sqlite database
            %   setting: either 'unshuffled', 'shuffled', or 'diagonal'
            %   sess: a string cell of session codes denoting which
            %       sessions to read out
            %   samp_size: how many samples to load, for a given set of
            %       parameters, e.g. 20 or 80
            bc = @DecodeTensor.build_command_sess;
            q = @Utils.cf_p;
            res = q(1,@(s)conn.fetch(bc(s, setting, 'MSE', [], 'max')), sess);
            mask = ~cellfun(@isempty, res);
            res = res(mask);
            n_sizes = q(1,@(r)double(cell2mat(r(:,1))), res);
            imse = q(1,@(r)1./cell2mat(r(:,3)), res);
            [n_sizes, imse] = Utils.cf_p2(1,...
                @(n,i)MultiSessionVisualizer.regroup(n, i, samp_size),...
                n_sizes, imse);
        end
        
        
        function plot_decoding_curve(sess, sp_, n_sizes, imse_s, I0_fit_s, N_fit_s, color, isrms)
            %plot a subset of the sessions as decoding curves + curve fits
            %Inputs:
            %   sess: cell array of strings of the session codes
            %   sp_: subset of the indices in sess to display
            %   n_sizes: x-axis values for each session in sess, as cell
            %   imse_s: y-axis values for each session in sess, as cell
            %   I0_fit_s: I0 fit value for each session in sess, as numeric array
            %   N_fit_s: N fit value for each session in sess, as numeric
            %       array
            %   color: The color of the curves, e.g. 'b', 'r', 'm'
            %   isrms (optional, default false): plot as RMS instead of
            %      IMSE, this applies @(x)x.^(-1/2) to the values in imse_s
            if ~exist('isrms', 'var')
                isrms = false;
            end
            if ~isrms
                pre = @(x)x;
                cap = 1;
            else
                pre = @(x)x.^(-1/2);
                cap = 0.4;
            end
            for j = 1:numel(sess)
                if ismember(j, sp_)
                    n = n_sizes{j};
                    %i = imse{i};
                    i_s = imse_s{j};
                    %plot(series_fits{2}{j}, 'r');
                    n_f = 1:500;
                    plot(n_f, pre(I0_fit_s(j).*n_f./(1 + n_f./N_fit_s(j))), color);
                    hold on;
                    errorbar(n, mean(pre(i_s)), std(pre(i_s))./sqrt(size(i_s,1)),...
                        color, 'Capsize', cap, 'LineStyle', 'none');
                end
            end
        end
        
        function adjr2 = plot_regress(x, y, x_conf, y_conf, mouse_names, color, varargin)
            p = inputParser;
            p.addOptional('xlim', [], @isnumeric);
            p.addOptional('text_coords', [], @isnumeric);
            p.addOptional('show_adjr2', true, @islogical);
            p.addOptional('fix_expo', false, @islogical);
            p.addParameter('dotsize', 4, @isnumeric);
            p.parse(varargin{:});
            
            %figure;
            dotsize = p.Results.dotsize;
            
            %InfoLimit = N_fit.*I0_fit;
            %InfoLimit_conf = abs(InfoLimit).*sqrt((N_conf./N_fit).^2 + (I0_conf./I0_fit).^2);
            
            errorbar(x, y, y_conf, y_conf,...
                x_conf, x_conf, 'LineStyle', 'none', 'Color', 'k', 'CapSize', 1);
            hold on;
            scatter(x, y, dotsize, DecodeTensor.mcolor(mouse_names, false), 'filled');
            
            xlim([min(x) - 0.1*range(x), max(x) + 0.1*range(x)]);
            ylim([min(y) - 0.1*range(y), max(y) + 0.1*range(y)]);
            %[pearson, pearson_p] = corr(x(:), y(:), 'type', 'Pearson');
            %[kendall, kendall_p] = corr(x(:), y(:), 'type', 'Kendall');
            %[spearman, spearman_p] = corr(x(:), y(:), 'type', 'Spearman');
            %fprintf('Pearson rho: %f, p=%f\nSpearman rho: %f, p=%f\nKendall tau: %f, p=%f\n',...
            %    pearson, pearson_p, spearman, spearman_p, kendall, kendall_p);
            [fitresult, adjr2] = Utils.regress_line(x, y);
            h_ = plot(fitresult); legend off
            h_.Color = color;
            h_.LineStyle = '--';
            if ~isempty(p.Results.xlim)
                xlim(p.Results.xlim);
            end
            if p.Results.show_adjr2
                if ~isempty(p.Results.text_coords)
                    text(p.Results.text_coords(1), p.Results.text_coords(2),...
                        sprintf('{\\itR}^2 = %.2f', adjr2));
                else
                    xl_ = xlim;
                    yl_ = ylim;
                    xl_l = [max(xl_(1),min(x)) min(xl_(2),max(x))];
                    yl_l = [max(yl_(1),min(y)) min(yl_(2),max(y))];
                    text(xl_l(1)+3/4*diff(xl_l), yl_l(1)+1/4*diff(yl_l),...
                        sprintf('{\\itR}^2 = %.2f', adjr2));
                end
            end
            
            if p.Results.fix_expo
                Utils.fix_exponent(gca, 'x', 0);
            end
            %xlabel 'IMSE limit'
            %ylabel '(d mu)^2 slope / sigma^2 slope'
        end
        
        function adjr2_list = plot_regress_series(x, y, x_conf, y_conf, mouse_list, color, varargin)
            p = inputParser;
            p.addOptional('xlim', [], @isnumeric);
            p.addOptional('text_coords', [], @isnumeric);
            p.addOptional('show_adjr2', true, @islogical);
            p.addOptional('fix_expo', false, @islogical);
            p.parse(varargin{:});
            
            mouse_names = unique(mouse_list);
            num_mice = numel(mouse_names);
            n_rows = round(sqrt(num_mice));
            n_cols = ceil(num_mice/n_rows);
            for m_i = 1:num_mice
                subplot(n_rows, n_cols, m_i);
                f_ = strcmp(mouse_names{m_i}, mouse_list);
                x_ = x(f_);
                y_ = y(f_);
                x_conf_ = x_conf(f_);
                y_conf_ = y_conf(f_);
                m_ = mouse_list(f_);
                
                adjr2_list(m_i) = PanelGenerator.plot_regress(x_, y_,...
                    x_conf_, y_conf_, m_, color, varargin{:});
                xlabel ''
                ylabel ''
                
                %xlim([0 500]);
                %ylim([0 max_y_val]);
                title(mouse_names{m_i});
            end %mice id
        end
        
        function adjr2 = plot_regress_averaged(x, y, x_conf, y_conf, mouse_list, color, varargin)
            p = inputParser;
            p.addOptional('xlim', [], @isnumeric);
            p.addOptional('text_coords', [], @isnumeric);
            p.addOptional('show_adjr2', true, @islogical);
            p.addParameter('dotsize', 4, @isnumeric);
            p.parse(varargin{:});
            
            mouse_names = unique(mouse_list);
            num_mice = numel(mouse_names);

            
            for m_i = 1:num_mice
                
                f_ = strcmp(mouse_names{m_i}, mouse_list);
                x_(m_i) = mean(x(f_));
                y_(m_i) = mean(y(f_));
                x_conf_(m_i) = sqrt(mean(x_conf(f_).^2));%sqrt(1.96.^2*var(x(f_)) + mean(x_conf(f_).^2));
                y_conf_(m_i) = sqrt(mean(y_conf(f_).^2));%sqrt(1.96.^2*var(y(f_)) + mean(y_conf(f_).^2));
                m_ = mouse_list(f_);
                
                %adjr2_list(m_i) = PanelGenerator.plot_regress(x_, y_,...
                %    x_conf_, y_conf_, m_, color, varargin{:});
                
                %%xlim([0 500]);
                %%ylim([0 max_y_val]);
                %title(mouse_names{m_i});
            end %mice id
            adjr2 = PanelGenerator.plot_regress(x_, y_, x_conf_, y_conf_, mouse_names, color, varargin{:});
        end
        
        function aux_decoding_curves(fname, sess, mouse_names, n_sizes, imse, imse_alt,...
                I0_fit, I0_fit_alt, N_fit, N_fit_alt, color, color_alt,...
                ylim_imse, ylim_rms, ylim_multi, non_inset)
            %For two series of IMSE curves (i.e. unshuffled-shuffled, or
            %unshuffled-diagonal), plot selected sessions as IMSE, RMS, along with the fit curve and
            %plot all sessions by mouse.
            %Inputs:
            %   fname: base filename with directory,
            %       'figure1_pdf/decoding_curves/fit_decoding_curves.pdf'
            %   sess: a char cell array of session IDs
            %   mouse_names: char cell array of mouse names corresponding
            %       to session IDs in sess
            %   n_sizes: x-axis values for each session in sess, as cell
            %   imse: y-axis values for each session in sess, as cell
            %   imse_alt: y-axis values for each session in sess, as cell,
            %       for the other condition
            %   I0_fit: I0 fit value for each session in sess, as numeric array
            %   I0_fit_alt: I0 fit value for each session in sess, as
            %       numeric array, for the other condition
            %   N_fit: N fit value for each session in sess, as numeric
            %       array
            %   N_fit_alt: N fit value for each session in sess, as numeric
            %       array, for the other condition
            %   color: The color of the curves, e.g. 'b', 'r', 'm' (should
            %       be 'b' for unshuffled)
            %   color_alt: The color of the curves, for the other condition
            %       e.g. 'b', 'r', 'm' (should be 'r' for shuffled, 'm' for diagonal)
            %   ylim_imse: The y axis limits for the IMSE curves
            %   ylim_rms: The y axis limits for the RMS curves
            %   ylim_multi: The common y axis limits for the complete set
            %       of curves. If given as [y1 y2] then y1 is ignored and
            %       assumed to be 0. Otherwise sets only the upper y limit.
            %       set to Inf for each mouse to have its own max
            figure('FileName', fname);
            show_mice = {'Mouse2022', 'Mouse2024', 'Mouse2028'};
            %[~,m_,sp_] = DecodeTensor.special_sess_id_list;<>
            %show_filter = ismember(m_, show_mice);
            %sp_ = sp_(show_filter);
            sp_ = SessManager.special_sessions_usable_index(show_mice);
            PanelGenerator.plot_decoding_curve(sess, sp_, n_sizes, imse_alt, I0_fit_alt, N_fit_alt, color_alt);
            hold on;
            PanelGenerator.plot_decoding_curve(sess, sp_, n_sizes, imse, I0_fit, N_fit, color);
            xlabel 'Number of cells'
            ylabel '1/MSE (cm^{-2})'
            %ylim([0 0.16]);
            ylim(ylim_imse);
            legend off
            figure_format([2 2.5]);
            Utils.printto;
            
            [pref, fn, ext] = fileparts(fname);
            figure('FileName', fullfile(pref, ['inset' fn ext]));
            PanelGenerator.plot_decoding_curve(sess, sp_, n_sizes, imse_alt, I0_fit_alt, N_fit_alt, color_alt, true);
            hold on;
            PanelGenerator.plot_decoding_curve(sess, sp_, n_sizes, imse, I0_fit, N_fit, color, true);
            l_ = refline(0, 5); l_.Color = 'k'; l_.LineStyle = ':';
            xlabel 'Number of cells'
            ylabel 'RMS Error (cm)'
            %ylim([2 50]);
            ylim(ylim_rms);
            set(gca, 'YScale', 'log');
            set(gca, 'YTick', [1 2 5 10 20 50]);
            legend off
            if exist('non_inset', 'var') && non_inset
                %figure_format;
                figure_format([2 2.5]);
            else
                figure_format('boxsize', [0.6 0.8], 'fontsize', 5);
            end
            Utils.printto;
            
            MultiSessionVisualizer.plot_series(n_sizes, {imse_alt, imse}, {color_alt, color}, mouse_names, ylim_multi);
            xlabel 'Number of cells'
            ylabel '1/MSE (cm^{-2})'
            multi_figure_format;
            Utils.printto('supplements_pdf/decoding_curves', ['multi_' fn ext]);
        end
        
        function aux_param_bns(fname, param, param_alt, param_conf, param_conf_alt,...
                mouse_names, x_labels, y_label, y_lim, fix_expo, log_scale, yticks)
            %Make both grouped ball-and-stick plots and a set of
            %ball-and-stick plots for each mouse comparing param and
            %param_alt in a paired comparison
            %Inputs:
            %   fname: base filename with directory, e.g.
            %       figure1_pdf/decoding_curves/grouped_I0_fit.pdf
            %   param: for each session, the value of a parameter fit, as
            %       numeric array
            %   param_alt: equivalent to param for the condition compared
            %       to
            %   param_conf: symmetric 95% confidence interval for param, as
            %       numeric array, given as "errorbar"
            %   param_conf_alt: equivalent to param_conf for the condition
            %       compared to
            %   mouse_names: char cell array of mouse names corresponding
            %       to session IDs in sess, and param
            %   x_labels: The names of the two conditions, as a cell
            %   y_label: The text for ylabel
            %   y_lim: y axis limits
            %   fix_expo: boolean, whether or not to fix exponents greater
            %       than 10^4, to make them valid for publication
            %   log_scale: boolean, log scale on the y axis or not
            [~, fn, ext] = fileparts(fname);
            figure('FileName', fname);
            Utils.bns_groupings(param, param_alt, param_conf, param_conf_alt, mouse_names, true, x_labels, log_scale);

            if ~isempty(y_lim)
                ylim(y_lim);
            end
            
            if exist('yticks', 'var') && ~isempty(yticks)
                set(gca, 'YTick', yticks);
            end
            
            ylabel(y_label);
            figure_format;
            if fix_expo
                Utils.fix_exponent(gca, 'Y', 0);
            end
            Utils.printto;
            
            figure;
            Utils.bns_groupings(param, param_alt, param_conf, param_conf_alt, mouse_names, false, x_labels, log_scale);
            
            if ~isempty(y_lim)
                ylim(y_lim);
            end
            
            ylabel(y_label);
            xlabel 'Mouse index';
            Utils.specific_format('MBNS');
            if fix_expo
                Utils.fix_exponent(gca, 'Y', 0);
            end
            Utils.printto('supplements_pdf/param_bns', ['multi_' fn ext]);
        end
        
        function aux_many_regress(x_cell, y_cell, x_conf_cell, y_conf_cell, mouse_list, x_lab_cell, y_lab_cell)
            n_rows = numel(x_cell);
            n_cols = numel(y_cell);
            assert(numel(x_conf_cell)==n_rows);assert(numel(y_conf_cell)==n_cols);assert(numel(x_lab_cell)==n_rows);assert(numel(y_lab_cell)==n_cols);
            figure;
            for r = 1:n_rows
                for c = 1:n_cols
                    subplot(n_rows, n_cols, sub2ind([n_cols n_rows], c,r));
                    my_x = x_cell{r};
                    my_x_conf = x_conf_cell{r};
                    my_x_lab = x_lab_cell{r};
                    my_y = y_cell{c};
                    my_y_conf = y_conf_cell{c};
                    my_y_lab = y_lab_cell{c};
                    
                    adjr2 = PanelGenerator.plot_regress(my_x, my_y, my_x_conf, my_y_conf, mouse_list, 'k');
                    xlabel(my_x_lab);
                    ylabel(my_y_lab);
                    title(sprintf('adj. R^2 = %.2f', adjr2));
                end
            end
            
        end
        
        function aux_regressions(fname, x, y, x_conf, y_conf, mouse_list,...
                color, text_coord, xlab, ylab, fix_expo, savedir_main, savedir_sup, special_xlim)
            [~, a_, b_] = fileparts(fname);
            fname = [a_ b_];
            
            figure;
            PanelGenerator.plot_regress(x, y,...
                x_conf, y_conf, mouse_list, color,...
                'text_coord', text_coord);
            xlabel(xlab);
            ylabel(ylab);
            if exist('special_xlim', 'var')
                xlim(special_xlim);
            end
            figure_format('factor', 1.6);
            if ~isequal(fix_expo, false)
                Utils.fix_exponent(gca, fix_expo, 0);
            end
            Utils.printto(savedir_main, fname);
            
            figure;
            PanelGenerator.plot_regress_averaged(x, y,...
                x_conf, y_conf, mouse_list, color,...
                'text_coord', text_coord);
            xlabel(xlab);
            ylabel(ylab);
            figure_format('factor', 1.6);
            if ~isequal(fix_expo, false)
                Utils.fix_exponent(gca, fix_expo, 0);
            end
            Utils.printto(savedir_sup, ['avg_' fname]);
            
            figure;
            PanelGenerator.plot_regress_series(x, y,...
                x_conf, y_conf, mouse_list, color);
            xlabel(xlab);
            ylabel(ylab);
            if fix_expo
                multi_figure_format('fix_expo', {'x', 0});
            else
                multi_figure_format;
            end
            Utils.printto(savedir_sup, ['multi_' fname]);
        end
        
        function val_reporter(v1, v2, c1, c2, lab1, lab2, mouse_names, post_func)
            ea = @(v,c) sqrt(1.96.^2.*var(v)/numel(v) + mean(c.^2));
            names_uniq = unique(mouse_names);
            
            [agg_v1, agg_v2, agg_c1, agg_c2] = deal(zeros(1, numel(names_uniq)));
            for m_i = 1:numel(names_uniq)
                m_name = names_uniq{m_i};
                filt = strcmp(mouse_names, m_name);
                
                agg_v1(m_i) = mean(v1(filt));
                agg_c1(m_i) = ea(v1(filt), c1(filt));
                
                agg_v2(m_i) = mean(v2(filt));
                agg_c2(m_i) = ea(v2(filt), c2(filt));
            end
            
            mv1 = mean(agg_v1);
            mv2 = mean(agg_v2);
            
            mc1 = ea(agg_v1, agg_c1);
            mc2 = ea(agg_v2, agg_c2);
            
            if ~exist('post_func', 'var')
                fprintf('Value of %s across mice: %e +- %e (95%% conf)\n', lab1, mv1, mc1);
                fprintf('Value of %s across mice: %e +- %e (95%% conf)\n', lab2, mv2, mc2);
            else
                mv1_avg = post_func(mv1);
                mv2_avg = post_func(mv2);
                mv1_lower = post_func(mv1 - mc1);
                mv2_lower = post_func(mv2 - mc2);
                mv1_upper = post_func(mv1 + mc1);
                mv2_upper = post_func(mv2 + mc2);
                fprintf('Value of %s, postfunc: %e - %e, mid: %e\n', lab1, mv1_lower, mv1_upper, mv1_avg);
                fprintf('Value of %s, postfunc: %e - %e, mid: %e\n', lab2, mv2_lower, mv2_upper, mv2_avg);
            end
        end
    end
    
    
    methods(Static)
        function decode_demo(varargin) %temporal_mixing
            p = inputParser;
            p.addParameter('temporal_mixing', true, @islogical);
            p.addParameter('ied', false, @islogical);
            p.addParameter('make_fig', true, @islogical);
            p.addParameter('pls', false, @islogical);
            p.parse(varargin{:});
            
            data_source = DecodeTensor.cons_filt(70, true); %was 70, 10, 34
            opt = DecodeTensor.default_opt;
            load(data_source{1});
            traces = tracesEvents.rawTraces;
            if p.Results.ied
                traces = iterative_event_detection(traces);
            end
            [~, ~, trial_start, trial_end, trial_dir, ~, ks] =...
                DecodeTensor.new_sel(tracesEvents.position, opt);
            within_trial = zeros(size(ks));
            trial_start = trial_start(trial_dir==1);
            trial_end = trial_end(trial_dir==1);
            n_trials = numel(trial_start);
            for i = 1:n_trials
                within_trial(trial_start(i):trial_end(i)) = i;
            end
            %%within_trial = within_trial .* (mod(ks,2)==1);
            %first_half = (within_trial <= n_trials/2) & (within_trial ~= 0);
            %second_half = (within_trial > n_trials/2) & (within_trial ~= 0);
            time_coord = (1:numel(ks))./opt.samp_freq;
            if p.Results.temporal_mixing
                first_half = (mod(within_trial,2) == 1) & (within_trial ~= 0);
                second_half = (mod(within_trial,2) == 0) & (within_trial ~= 0);
            else
                first_half = (within_trial <= n_trials/2) & (within_trial ~= 0);
                second_half = (within_trial > n_trials/2) & (within_trial ~= 0);
            end
            
            ks_first_half = ks(first_half);
            ks_second_half = ks(second_half);
            X_first_half = traces(first_half, :);
            X_second_half = traces(second_half, :);
            X_first_half_s = shuffle(X_first_half, ks_first_half);
            X_second_half_s = shuffle(X_second_half, ks_second_half);
            
            alg = my_algs('ecoclin');
            model1 = alg.train(X_first_half, ks_first_half); disp('m1');
            model1_s = alg.train(X_first_half_s, ks_first_half); disp('m1s');
            model2 = alg.train(X_second_half, ks_second_half); disp('m2');
            model2_s = alg.train(X_second_half_s, ks_second_half); disp('m2s');
            
            ps_first_half = alg.test(model2, X_first_half);
            ps_first_half_s = alg.test(model2_s, X_first_half_s);
            ps_first_half_d = alg.test(model2_s, X_first_half);
            ps_second_half = alg.test(model1, X_second_half);
            ps_second_half_s = alg.test(model1_s, X_second_half_s);
            ps_second_half_d = alg.test(model1_s, X_second_half);
            
            me_ = @(k,p) mean(abs(ceil(k(:)/2)-ceil(p(:)/2))).*opt.bin_width;
            mse_ = @(k,p) mean((ceil(k(:)/2)-ceil(p(:)/2)).^2).*opt.bin_width.^2;
            fprintf('Mean error: unsh: %f, sh: %f\n',  me_(ks_first_half, ps_first_half),...
                me_(ks_first_half, ps_first_half_s));
            fprintf('Mean error: unsh: %f, sh: %f\n',  me_(ks_second_half, ps_second_half),...
                me_(ks_second_half, ps_second_half_s));
            
            
            tr_track = within_trial(first_half|second_half);
            %ks_ = [ks_first_half; ks_second_half];
            %ps_ = [ps_first_half; ps_second_half];
            %ps_s_ = [ps_first_half_s; ps_second_half_s];
            ks_ = ks(within_trial~=0);
            ps_(first_half) = ps_first_half;
            ps_(second_half) = ps_second_half;
            ps_ = ps_(within_trial~=0)';
            ps_s_(first_half) = ps_first_half_s;
            ps_s_(second_half) = ps_second_half_s;
            ps_s_ = ps_s_(within_trial~=0)';
            
            for i = 1:max(tr_track)
                tr_err(i) = me_(ks_(tr_track==i), ps_(tr_track==i));
                tr_err_s(i) = me_(ks_(tr_track==i), ps_s_(tr_track==i));
            end
            
            
            figure('FileName', 'figure1_pdf\demo\decoding_demo.pdf');
            t = (1:numel(ks_first_half))./opt.samp_freq;
            t_start = 83.5 + 2.05; %which to show
            t_end = 83.5 + 4.05;%90;
            
            hold on;
            h(1) = plot(t - t_start, (ceil(ps_first_half/2) - 0.5)*opt.bin_width, '-b');
            h(2) = plot(t - t_start, (ceil(ps_first_half_s/2) - 0.5)*opt.bin_width, '-r');
            h(3) = plot(t - t_start, (ceil(ks_first_half/2) - 0.5)*opt.bin_width, '-k');
            
            trial_boundaries = (find(diff(within_trial(first_half))>0)+0.5)./opt.samp_freq - t_start;
            for i = 1:numel(trial_boundaries)
                x_ = trial_boundaries(i);
                line([x_ x_], ylim, 'Color', 'k', 'LineStyle', '--');
            end
            xlim([t_start t_end] - t_start);
            ylim([-Inf Inf]);
            xlabel 'Time (s)';
            ylabel 'Position (cm)';
            legend(h, 'Real', 'Shuffled', 'Place bin');
            legend boxoff
            figure_format('boxsize', [0.8 0.7]*1.05); box on;
            if p.Results.make_fig
                Utils.printto;
            end
            
            figure('FileName', 'figure1_pdf\demo\decoding_demo_diagonal.pdf');
            t = (1:numel(ks_first_half))./opt.samp_freq;
            %t_start = 83.5;
            %t_end = 90;
            
            hold on;
            h(1) = plot(t - t_start, (ceil(ps_first_half/2) - 0.5)*opt.bin_width, '-b');
            h(2) = plot(t - t_start, (ceil(ps_first_half_d/2) - 0.5)*opt.bin_width, '-m');
            h(3) = plot(t - t_start, (ceil(ks_first_half/2) - 0.5)*opt.bin_width, '-k');
            
            trial_boundaries = (find(diff(within_trial(first_half))>0)+0.5)./opt.samp_freq - t_start;
            for i = 1:numel(trial_boundaries)
                x_ = trial_boundaries(i);
                line([x_ x_], ylim, 'Color', 'k', 'LineStyle', '--');
            end
            xlim([t_start t_end] - t_start);
            ylim([-Inf Inf]);
            xlabel 'Time (s)';
            ylabel 'Position (cm)';
            legend(h, 'Real', 'Diagonal', 'Place bin');
            legend boxoff
            figure_format('boxsize', [0.8 0.7]*1.05); box on;
            if p.Results.make_fig
                Utils.printto;
            end
            
            if p.Results.pls
                [~, stats] = Utils.pls_plot([X_first_half;X_second_half],...
                    [time_coord(first_half)', ceil(ks_first_half/2), mod(ks_first_half,2);...
                    time_coord(second_half)', ceil(ks_second_half/2),mod(ks_second_half,2)]);
                Utils.pls_plot([X_first_half_s;X_second_half_s],...
                    [time_coord(first_half)', ceil(ks_first_half/2), mod(ks_first_half,2);...
                    time_coord(second_half)', ceil(ks_second_half/2),mod(ks_second_half,2)], 'stats', stats, 'xl_', xlim, 'yl_', ylim);
            end
        end
        
        function confusion(remake)
            if ~exist('remake', 'var')
                remake = false;
            end
            savedir = 'figure1_pdf/confusion';
            if ~exist(savedir, 'dir')
                mkdir(savedir);
            end
            
            ap = @(x) fullfile(savedir, x);
            
            
            load('confusion_single_session_agg_190730-145852_0.mat');
            nt = res(1).num_trials;
            C_pct = cat(3, res.C)/nt*100;
            C_s_pct = cat(3, res.C_s)/nt*100;
            C_d_pct = cat(3, res.C_d)/nt*100;
            
            
            fname = ap('confusion_unshuffled.pdf');
            if remake || ~exist(fname, 'file')
                figure('FileName', fname);
                PanelGenerator.plot_confusion(C_pct, 'Confusion (%)');
                colormap gray
                colormap(flipud(colormap(gca)));
                Utils.printto;
            end
            
            fname = ap('confusion_shuffled.pdf');
            if remake || ~exist(fname, 'file')
                figure('FileName', fname);
                PanelGenerator.plot_confusion(C_s_pct, 'Confusion (%)');
                colormap gray
                colormap(flipud(colormap(gca)));
                Utils.printto;
            end
            
            fname = ap('confusion_diagonal.pdf');
            if remake || ~exist(fname, 'file')
                disp('making diagonal');
                figure('FileName', fname);
                PanelGenerator.plot_confusion(C_d_pct, 'Confusion (%)');
                colormap gray
                colormap(flipud(colormap(gca)));
                Utils.printto;
            end
            
            fname = ap('confusion_diff.pdf');
            if remake || ~exist(fname, 'file')
                figure('FileName', fname);
                PanelGenerator.plot_confusion(C_s_pct - C_pct, 'Confusion difference (%)', false);
                colormap(bluewhitered);
                Utils.printto;
            end
            
            fname = ap('bin_error_rate.pdf');
            if remake || ~exist(fname, 'file')
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
            
            fname = ap('bin_error_rate_diagonal.pdf');
            if remake || ~exist(fname, 'file')
                figure('FileName', fname);
                m = 100-diag(squeeze(mean(C_pct,3)));
                e = diag(squeeze(std(C_pct,[],3)))./sqrt(size(C_pct,3));
                errorbar(m, e, 'b', 'CapSize', 1);
                xlabel 'Place bin'
                ylabel 'Error rate (%)'
                hold on
                m_d = 100-diag(squeeze(mean(C_d_pct,3)));
                e_d = diag(squeeze(std(C_d_pct,[],3)))./sqrt(size(C_d_pct,3));
                errorbar(m_d, e_d, 'm', 'CapSize', 1);
                errorbar(m_d - m, sqrt(e.^2 + e_d.^2), 'k', 'DisplayName', 'Difference', 'CapSize', 1);
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
                samp_size = 20;
                %[sess, mouse_names] = DecodeTensor.filt_sess_id_list;
                [sess, mouse_names] = SessManager.usable_sess_id_list;
                [n_sizes, imse] = PanelGenerator.db_imse_reader(conn, 'unshuffled', sess, samp_size);
                [n_sizes_s, imse_s] = PanelGenerator.db_imse_reader(conn, 'shuffled', sess, samp_size);
                [n_sizes_d, imse_d] = PanelGenerator.db_imse_reader(conn, 'diagonal', sess, samp_size);
                assert(isequal(n_sizes, n_sizes_s), 'mismatch between unshuffled and shuffled sampling');
                assert(isequal(n_sizes, n_sizes_d), 'mismatch between unshuffled and diagonal sampling');
                
                [series_fits{1}, series_gof{1}] = Utils.cf_p2(2,@(n,m)createFit_infoSaturation(n(:),mean(m)'), n_sizes, imse);
                progressbar(1/3/2);
                [series_fits{2}, series_gof{2}] = Utils.cf_p2(2,@(n,m)createFit_infoSaturation(n(:),mean(m)'), n_sizes, imse_s);
                progressbar(2/3/2);
                [series_fits{3}, series_gof{3}] = Utils.cf_p2(2,@(n,m)createFit_infoSaturation(n(:),mean(m)'), n_sizes, imse_d);
                progressbar(3/3/2);
                [series_exp_fits{1}, series_exp_gof{1}] = Utils.cf_p2(2,@(n,m)createFit_exp(n(:),mean(m)'), n_sizes, imse);
                progressbar(4/3/2);
                [series_exp_fits{2}, series_exp_gof{2}] = Utils.cf_p2(2,@(n,m)createFit_exp(n(:),mean(m)'), n_sizes, imse_s);
                progressbar(5/3/2);
                [series_exp_fits{3}, series_exp_gof{3}] = Utils.cf_p2(2,@(n,m)createFit_exp(n(:),mean(m)'), n_sizes, imse_d);
                progressbar(6/3/2);
                
                [I0_fit, I0_conf] = Utils.fit_get(series_fits{1}, 'I_0');
                [I0_fit_s, I0_conf_s] = Utils.fit_get(series_fits{2}, 'I_0');
                [I0_fit_d, I0_conf_d] = Utils.fit_get(series_fits{3}, 'I_0');
                
                [N_fit, N_conf] = Utils.fit_get(series_fits{1}, 'N');
                [N_fit_s, N_conf_s] = Utils.fit_get(series_fits{2}, 'N');
                [N_fit_d, N_conf_d] = Utils.fit_get(series_fits{3}, 'N');
                
                save(save_file, 'sess', 'mouse_names', 'n_sizes',...
                    'imse', 'imse_s', 'imse_d', 'series_fits', 'series_gof', 'I0_fit', 'I0_conf',...
                    'I0_fit_s', 'I0_conf_s', 'N_fit', 'N_conf',...
                    'N_fit_s', 'N_conf_s', 'series_exp_fits', 'series_exp_gof',...
                    'I0_fit_d', 'I0_conf_d', 'N_fit_d', 'N_conf_d');
            else
                load(save_file);
            end
            
            r2 = cellfun(@(x)x.rsquare, series_gof{1});
            r2_s = cellfun(@(x)x.rsquare, series_gof{2});
            r2_d = cellfun(@(x)x.rsquare, series_gof{3});
            
            r2_exp = cellfun(@(x)x.rsquare, series_exp_gof{1});
            r2_s_exp = cellfun(@(x)x.rsquare, series_exp_gof{2});
            r2_d_exp = cellfun(@(x)x.rsquare, series_exp_gof{3});
            
            fprintf('Unshuf R^2: %f - %f, median: %f\n', min(r2), max(r2), median(r2));
            fprintf('Shuf R^2: %f - %f, median: %f\n', min(r2_s), max(r2_s), median(r2_s));
            fprintf('Diag R^2: %f - %f, median: %f\n', min(r2_d), max(r2_d), median(r2_d));
            fprintf('Unshuf (exp) R^2: %f - %f, median: %f\n', min(r2_exp), max(r2_exp), median(r2_exp));
            fprintf('Shuf (exp) R^2: %f - %f, median: %f\n', min(r2_s_exp), max(r2_s_exp), median(r2_s_exp));
            fprintf('Diag (exp) R^2: %f - %f, median: %f\n', min(r2_d_exp), max(r2_d_exp), median(r2_d_exp));
            
            
            fname = ap('decoding_curve_fit.pdf');
            if p.Results.remake || ~exist(fname, 'file')
                PanelGenerator.aux_decoding_curves(fname, sess, mouse_names, n_sizes, imse, imse_s,...
                    I0_fit, I0_fit_s, N_fit, N_fit_s, 'b', 'r',...
                    [0 0.16], [2 50], [0 Inf], true);
            end
            
            fname = ap('decoding_curve_fit_diagonal');
            if p.Results.remake || ~exist(fname, 'file')
                PanelGenerator.aux_decoding_curves(fname, sess, mouse_names, n_sizes, imse, imse_d,...
                    I0_fit, I0_fit_d, N_fit, N_fit_d, 'b', 'm',...
                    [0 0.05], [4 50], [0 Inf], true);
            end
            
            fname = ap('grouped_I0_fit.pdf');
            if p.Results.remake || ~exist(fname, 'file')
                %figure('FileName', fname);
                %Utils.bns_groupings(I0_fit, I0_fit_s, I0_conf, I0_conf_s, mouse_names, true);
                %ylim([-Inf Inf]);
                %ylabel(sprintf('I_0 fit value\n(cm^{-2}neuron^{-1})'));
                %figure_format;
                %Utils.fix_exponent(gca, 'Y', 0);
                %Utils.printto;
                PanelGenerator.aux_param_bns(fname, I0_fit, I0_fit_s, I0_conf, I0_conf_s, mouse_names,...
                    {'Unshuffled', 'Shuffled'}, sprintf('I_0 fit value\n(cm^{-2}neuron^{-1})'), [-Inf Inf], true, false);
                PanelGenerator.val_reporter(I0_fit, I0_fit_s, I0_conf, I0_conf_s, 'I_0', 'I_0 (shuf)', mouse_names);
                
            end
            
            
            fname = ap('grouped_I0N_fit.pdf');
            if p.Results.remake || ~exist(fname, 'file')
                InfoLimit = N_fit.*I0_fit;
                InfoLimit_conf = abs(InfoLimit).*sqrt((N_conf./N_fit).^2 + (I0_conf./I0_fit).^2);
                
                InfoLimit_s = N_fit_s.*I0_fit_s;
                InfoLimit_conf_s = abs(InfoLimit_s).*sqrt((N_conf_s./N_fit_s).^2 + (I0_conf_s./I0_fit_s).^2);
                
                PanelGenerator.aux_param_bns(fname, InfoLimit, InfoLimit_s, InfoLimit_conf, InfoLimit_conf_s, mouse_names,...
                    {'Unshuffled', 'Shuffled'}, sprintf('I_0N fit value\n(cm^{-2})'), [1e-4 Inf], false, true, [1e-4 1e-2 1e0 1e2]);
                PanelGenerator.val_reporter(InfoLimit, InfoLimit_s, InfoLimit_conf, InfoLimit_conf_s, 'I_0N', 'I_0N (shuf)', mouse_names);
                PanelGenerator.val_reporter(log(InfoLimit), log(InfoLimit_s), InfoLimit_conf./InfoLimit, InfoLimit_conf_s./InfoLimit_s, 'log(I_0N)', 'log(I_0N) (shuf)', mouse_names, @exp);
                
            end
            
            fname = ap('grouped_I0_fit_diag.pdf');
            if p.Results.remake || ~exist(fname, 'file')
                %figure('FileName', fname);
                %Utils.bns_groupings(I0_fit, I0_fit_d, I0_conf, I0_conf_d, mouse_names, true, {'Unshuffled', 'Diagonal'});
                %ylim([-Inf Inf]);
                %ylabel(sprintf('I_0 fit value\n(cm^{-2}neuron^{-1})'));
                %figure_format;
                %Utils.fix_exponent(gca, 'Y', 0);
                %Utils.printto;
                PanelGenerator.aux_param_bns(fname, I0_fit, I0_fit_d, I0_conf, I0_conf_d, mouse_names,...
                    {'Unshuffled', 'Diagonal'}, sprintf('I_0 fit value\n(cm^{-2}neuron^{-1})'), [-Inf Inf], true, false);
                PanelGenerator.val_reporter(I0_fit, I0_fit_d, I0_conf, I0_conf_d, 'I_0', 'I_0 (diag)', mouse_names);
            end
            
            
            fname = ap('grouped_N_fit.pdf');
            if p.Results.remake || ~exist(fname, 'file')
                %figure('FileName', fname);
                %Utils.bns_groupings(N_fit, N_fit_s, N_conf, N_conf_s, mouse_names, true);
                %set(gca, 'YScale', 'log');
                %%ylim([-Inf Inf]);
                %ylabel(sprintf('N fit value\n(neuron)'));
                %figure_format;
                %Utils.printto;
                PanelGenerator.aux_param_bns(fname, N_fit, N_fit_s, N_conf, N_conf_s, mouse_names,...
                    {'Unshuffled', 'Shuffled'}, sprintf('N fit value\n(neuron)'), [], false, true);
                PanelGenerator.val_reporter(N_fit, N_fit_s, N_conf, N_conf_s, 'N', 'N (shuf)', mouse_names);
                PanelGenerator.val_reporter(log(N_fit), log(N_fit_s), N_conf./N_fit, N_conf_s./N_fit_s, 'log(N)', 'log(N) (shuf)', mouse_names, @exp);
            end
            
            fname = ap('grouped_n99_fit.pdf');
            if p.Results.remake || ~exist(fname, 'file')
                PanelGenerator.aux_param_bns(fname, 99*N_fit, 99*N_fit_s, 99*N_conf, 99*N_conf_s, mouse_names,...
                    {'Unshuffled', 'Shuffled'}, sprintf('Size at 99%% max.\n(neuron)'), [100 Inf], false, true, [1e2 1e4 1e6 1e8]);
                PanelGenerator.val_reporter(99*N_fit, 99*N_fit_s, 99*N_conf, 99*N_conf_s, '99N', '99N (shuf)', mouse_names);
                PanelGenerator.val_reporter(log(99*N_fit), log(99*N_fit_s), N_conf./N_fit, N_conf_s./N_fit_s, 'log(99N)', 'log(99N) (shuf)', mouse_names, @exp);
            end
            
            fname = ap('grouped_N_fit_diag.pdf');
            if p.Results.remake || ~exist(fname, 'file')
                %figure('FileName', fname);
                %Utils.bns_groupings(N_fit, N_fit_d, N_conf, N_conf_d, mouse_names, true, {'Unshuffled', 'Diagonal'});
                %set(gca, 'YScale', 'log');
                %ylabel(sprintf('N fit value\n(neuron)'));
                %figure_format;
                %Utils.printto;
                PanelGenerator.aux_param_bns(fname, N_fit, N_fit_d, N_conf, N_conf_d, mouse_names,...
                    {'Unshuffled', 'Diagonal'}, sprintf('N fit value\n(neuron)'), [], false, true);
                PanelGenerator.val_reporter(N_fit, N_fit_d, N_conf, N_conf_d, 'N', 'N (diag)', mouse_names);
                PanelGenerator.val_reporter(log(N_fit), log(N_fit_d), N_conf./N_fit, N_conf_d./N_fit_d, 'log(N)', 'log(N) (diag)', mouse_names, @exp);
            end
        end
        
        function medload
            load('MedLoad_agg_190705-171806_0.mat');
            n_sizes = {res.n_sizes};
            series = {{res.median_loadings}, {res.median_loadings_s}};
            
            series = Utils.cf_(@(m)Utils.cf_(@(x)max(x,[],3),m), series);
            mouse_name = {res.mouse_name};
            
            show_mice = {'Mouse2022', 'Mouse2024', 'Mouse2028'};
            
            %[~,m_,sp_] = DecodeTensor.special_sess_id_list;<>
            %show_filter = ismember(m_, show_mice);
            %sp_ = sp_(show_filter);
            %m_ = m_(show_filter);
            sp_ = SessManager.special_sessions_usable_index(show_mice);
            m_ = show_mice;
            figure('FileName', 'supplements_pdf/medload/medload_rasters.pdf');
            colorscale = 'log';
            for i = 1:numel(sp_)
                subplot(1,numel(sp_)+1, i);
                mean_median_loadings = squeeze(mean(abs(res(sp_(i)).median_loadings)));
                min_d = 30;
                ns = res(sp_(i)).n_sizes;
                im_data = (mean_median_loadings(ns >= min_d,1:min_d));
                %padded_im_data = nan(16 - size(im_data,1), size(im_data,2));
                %imagesc([im_data;padded_im_data], [-1.6 log10(0.3)]);
                surf(1:min_d, ns(ns>=min_d), im_data, 'EdgeColor', 'none');
                view(2);
                
                %set(gca, 'XScale', 'log');
                %set(gca, 'YScale', 'log');
                set(gca, 'ColorScale', colorscale);
                caxis([0.03 0.25]);
                xlim([1 min_d]);
                ylim([min_d+10, 500]);
                xlabel 'Fluctuation mode, i'
                ylabel 'Number of cells'
                title(sprintf('Mouse %s', m_{i}(end-1:end)), 'FontName', 'Helvetica', 'FontSize', 6, 'FontWeight', 'normal', 'Color', 'b');
                
                set(gca, 'FontSize', 6);
                set(gca, 'FontName', 'Helvetica');
                set(gca, 'TickLength', [0.02 0.02]);
                colorbar;
                %set(gca, 'YTick', [1 2 4 8 16 32 64].*min_d);
                %rectangle('Position',...
                %    0.5+[0 size(im_data,1) size(im_data,2) (48 - size(im_data,1))],...
                %    'FaceColor', 'w', 'EdgeColor', 'k', 'LineStyle', 'none');
                %set(gca, 'YTickLabel', 10*cellfun(@str2num, get(gca, 'YTickLabel')));
                box off;
                if i > 1
                    box off
                    xlabel ''
                    ylabel ''
                    set(gca, 'YTick', []);
                end
                %colorbar;
            end
            subplot(1, numel(sp_)+1, numel(sp_)+1);
            mean_median_loadings_s = squeeze(mean(abs(res(sp_(1)).median_loadings_s)));
            min_d = 30;
            ns = res(sp_(1)).n_sizes;
            im_data = (mean_median_loadings_s(ns >= min_d,1:min_d));
            surf(1:min_d, ns(ns>=min_d), im_data, 'EdgeColor', 'none');
            view(2);
            %set(gca, 'YScale', 'log');
            set(gca, 'ColorScale', colorscale);
            caxis([0.03 0.25]);
            xlim([1 min_d]);
            ylim([min_d+10, 500]);
            %xlabel 'Fluctuation mode, i'
            %ylabel 'Number of cells'
            title('Shuffled', 'FontName', 'Helvetica', 'FontSize', 6, 'FontWeight', 'normal', 'Color', 'r');
            
            set(gca, 'FontSize', 6);
            set(gca, 'FontName', 'Helvetica');
            set(gca, 'TickLength', [0.02 0.02]);
            set(gca, 'YTick', [1 2 4 8 16 32 64].*min_d);
            box off
            %axis off
            xlabel ''; ylabel '';
            set(gca, 'YTick', []);
            colorbar;
            set(gcf, 'Units', 'inches');
            set(gcf, 'Position', [8.5521    6.2292    8.3125    1.6146]);
            colormap parula;
            Utils.printto;
            
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
            
            MultiSessionVisualizer.plot_series(n_sizes, series, {'b','r'}, mouse_name);
            axs = findall(gcf, 'type', 'axes');
            set(axs, 'YScale', 'log');
            set(axs, 'XScale', 'log');
            xlabel 'Number of cells'
            ylabel 'max_i|cos(PC_i, Dm)|'
            multi_figure_format;
            Utils.printto('supplements_pdf/medload', 'multi_medload_curves.pdf');
            
            n_c = 50;
            fit_func = @(x,y)fit(log10(x(x>=n_c))',log10(mean(y(:,x>=n_c)))', 'poly1');
            [fr_, gf_] = cellfun(fit_func, n_sizes, series{1}, 'UniformOutput', false);
            [fr_s, gf_s] = cellfun(fit_func, n_sizes, series{2}, 'UniformOutput', false);
            
            gf_ = cell2mat(gf_);
            rsquare = [gf_.rsquare];
            fprintf('Unshuf: range %f-%f, median %f\n', min(rsquare), max(rsquare), median(rsquare));
            
            gf_s = cell2mat(gf_s);
            rsquare_s = [gf_s.rsquare];
            fprintf('Shuf: range %f-%f, median %f\n', min(rsquare_s), max(rsquare_s), median(rsquare_s));
            
            [rate_f, rate_f_conf] = Utils.fit_get(fr_, 'p1');
            [rate_f_s, rate_f_s_conf] = Utils.fit_get(fr_s, 'p1');
            figure('FileName', 'figure2_pdf/medload/inset.pdf');
            Utils.bns_groupings(rate_f, rate_f_s, rate_f_conf, rate_f_s_conf, mouse_name, true);
            hold on;
            %line(xlim-0.5, [0 0], 'Color', 'k', 'LineStyle', '-');
            line(xlim, [-0.5 -0.5], 'Color', 'k', 'LineStyle', ':');
            ylabel 'Fit exponent'
            ylim([-Inf Inf]);
            set(gca, 'XTickLabels', {'Real', 'Shuf.'});
            %figure_format([0.8 1]/2, 'fontsize', 4);
            Utils.specific_format('inset');
            Utils.printto;
            
            figure('FileName', 'supplements_pdf/medload/multi_inset.pdf');
            Utils.bns_groupings(rate_f, rate_f_s, rate_f_conf, rate_f_s_conf, mouse_name, false, {'Unshuffled', 'Shuffled'});
            hold on;
            line(xlim-0.5, [0 0], 'Color', 'k', 'LineStyle', '-');
            line(xlim, [-0.5 -0.5], 'Color', 'k', 'LineStyle', ':');
            ylabel 'Fit exponent'
            ylim([-Inf Inf]);
            %set(gca, 'XTickLabels', {'Unsh.', 'Sh.'});
            %figure_format([0.8 1]/2, 'fontsize', 4);
            Utils.specific_format('MBNS');
            Utils.printto;
            
            %PABLO POINT 4
            if false
                fprintf('The following are the values (there is no variation for max cells)\n for max_i cos(PC_i, Dm), at the maximal # of cells\n');
                max_overlap = cellfun(@(x)mean(x(:,end)), series{1}); %vector of length 107
                disp(max_overlap);
                
                if false
                dbfile = 'decoding_all_sess.db';
                conn = sqlite(dbfile);
                samp_size = 80;
                [sess, mouse_names] = DecodeTensor.filt_sess_id_list;
                [n_sizes, imse] = PanelGenerator.db_imse_reader(conn, 'unshuffled', sess, samp_size);
                [n_sizes_s, imse_s] = PanelGenerator.db_imse_reader(conn, 'shuffled', sess, samp_size);
                [n_sizes_d, imse_d] = PanelGenerator.db_imse_reader(conn, 'diagonal', sess, samp_size);
                assert(isequal(n_sizes, n_sizes_s), 'mismatch between unshuffled and shuffled sampling');
                assert(isequal(n_sizes, n_sizes_d), 'mismatch between unshuffled and diagonal sampling');
                
                max_imse = cell2mat(cellfun(@(x)x(:,end),imse,'UniformOutput', false));
                max_imse_s = cell2mat(cellfun(@(x)x(:,end),imse_s,'UniformOutput', false));
                
                mean_max_imse = mean(max_imse);
                mean_max_imse_s = mean(max_imse_s);
                
                sem_max_imse = sem(max_imse);
                sem_max_imse_s = sem(max_imse_s);
                
                max_num_neurons = cellfun(@(x)x(end),n_sizes);
                
                victors_suggestion = (mean_max_imse_s - mean_max_imse)./max_num_neurons;
                victors_suggestion_sem = sqrt( (sem_max_imse_s.^2 + sem_max_imse.^2)./max_num_neurons );
                
                
                [~, mouse_names] = DecodeTensor.filt_sess_id_list;
                figure;
                PanelGenerator.plot_regress(max_overlap, victors_suggestion,...
                    0*max_overlap, victors_suggestion_sem, mouse_names, 'g', 'text_coords', [0.25 -1e-4]);
                xlabel 'max_i|cos(PC_i, \Delta\mu)|'
                ylabel(sprintf('(IMSE_{sh} - IMSE_{re}) / (max n)\n(cm^{-2}%sneuron^{-1})',Utils.dot));
                Utils.fix_exponent(gca, 'y', 0);
                figure_format('factor', 1.6);
                Utils.printto('supplements_pdf', 'Victors_suggestion.pdf');
                
                figure;
                PanelGenerator.plot_regress_averaged(max_overlap, victors_suggestion,...
                    0*max_overlap, victors_suggestion_sem, mouse_names, 'g', 'text_coords', [0.25 -1e-4]);
                xlabel 'max_i|cos(PC_i, \Delta\mu)|'
                ylabel(sprintf('(IMSE_{sh} - IMSE_{re}) / (max n)\n(cm^{-2}%sneuron^{-1})',Utils.dot));
                Utils.fix_exponent(gca, 'y', 0);
                figure_format('factor', 1.6);
                Utils.printto('supplements_pdf', 'Victors_suggestion_grouped.pdf');
                
                pablos_suggestion_x = 1./sqrt(max_num_neurons) - max_overlap;
                pablos_suggestion_y = mean_max_imse_s - mean_max_imse;
                pablos_suggestion_y_sem = sqrt(sem_max_imse_s.^2 + sem_max_imse.^2);
                
                [~, mouse_names] = DecodeTensor.filt_sess_id_list;
                figure;
                PanelGenerator.plot_regress(pablos_suggestion_x, pablos_suggestion_y,...
                    0*pablos_suggestion_x, pablos_suggestion_y_sem, mouse_names, 'g', 'text_coords', [-0.22 0.1]);
                xlabel '1/sqrt(n) - max_i|cos(PC_i, \Delta\mu)|'
                ylabel(sprintf('(IMSE_{sh} - IMSE_{re}) (cm^{-2})'));
                %Utils.fix_exponent(gca, 'y', 0);
                figure_format('factor', 1.6);
                Utils.printto('supplements_pdf', 'Pablos_suggestion.pdf');
                
                figure;
                PanelGenerator.plot_regress_averaged(pablos_suggestion_x, pablos_suggestion_y,...
                    0*pablos_suggestion_x, pablos_suggestion_y_sem, mouse_names, 'g', 'text_coords', [-0.27 0.1]);
                xlabel '1/sqrt(n) - max_i|cos(PC_i, \Delta\mu)|'
                ylabel(sprintf('(IMSE_{sh} - IMSE_{re}) (cm^{-2})'));
                %Utils.fix_exponent(gca, 'y', 0);
                figure_format('factor', 1.6);
                Utils.printto('supplements_pdf', 'Pablos_suggestion_grouped.pdf');
                end
                
                
                if true
                    fit_savefile = 'decoding_curves_fits.mat';
                    recompute = false;
                    if recompute || ~exist(fit_savefile, 'file')
                        PanelGenerator.decoding_curves('remake', true, 'recompute', recompute);
                    end
                    load(fit_savefile);

                    if false
                    fprintf('The following is Pablo''s formula: (I0N(shuf)-I0N(real))\n');
                    pablos_formula = (I0_fit_s.*N_fit_s - I0_fit.*N_fit);
                    disp(pablos_formula);
                    figure;
                    scatter(max_overlap, pablos_formula);
                    xlabel 'max_i |cos(PC_i, \Delta\mu)|'
                    ylabel 'I_0N_{sh}-I_0N_{re}'
                    refline
                    [res_, GOF_] = fit(max_overlap(:), pablos_formula(:), 'poly1');
                    text(0.25, -0.6, sprintf('\\it{R}^2 = %.2e', GOF_.rsquare));

                    [~, mouse_names] = DecodeTensor.filt_sess_id_list;
                    figure;
                    PanelGenerator.plot_regress_averaged(max_overlap, pablos_formula,...
                        0*max_overlap, 0*pablos_formula, mouse_names, 'g',...
                        'xlim', [0.08 0.37], 'text_coords', [0.25 -5]);
                    xlabel 'max_i |cos(PC_i, \Delta\mu)|'
                    ylabel 'I_0N_{sh}-I_0N_{re}'
                    figure_format('factor', 1.6);
                    Utils.printto('supplements_pdf', 'Pablo_point_4.pdf')
                    end
                    
                    [~, mouse_names] = DecodeTensor.filt_sess_id_list;
                    figure;
                    PanelGenerator.plot_regress(max_overlap, N_fit,...
                        0*max_overlap, N_conf, mouse_names, 'g', 'show_adjr2', false);
                    xlabel 'max_i|cos(PC_i, \Delta\mu)|'
                    ylabel('N fit value');
                    %Utils.fix_exponent(gca, 'y', 0);
                    set(gca, 'YScale', 'log');
                    figure_format('factor', 1.6);
                    Utils.printto('supplements_pdf', 'Connection_to_N.pdf');
                    
                    %figure;
                    %PanelGenerator.plot_regress_averaged(max_overlap, N_fit,...
                    %    0*max_overlap, N_conf, mouse_names, 'g', 'text_coords', [0.22 100]);
                    %xlabel 'max_i|cos(PC_i, \Delta\mu)|'
                    %ylabel('N fit value');
                    %%Utils.fix_exponent(gca, 'y', 0);
                    %figure_format('factor', 1.6);
                    %Utils.printto('supplements_pdf', 'Connection_to_N_grouped.pdf');
                end
            end
        end
        
        
        function medload_with_mean
            load('MedLoad_agg_190705-171806_0.mat');
            n_sizes = {res.n_sizes};
            series = {{res.median_loadings}, {res.median_loadings_s}};
            

            series = Utils.cf_(@(m)Utils.cf_(@PanelGenerator.mean_func_trunc5,m), series); %using mean rather than max
            mouse_name = {res.mouse_name};
            %{
            show_mice = {'Mouse2022', 'Mouse2024', 'Mouse2028'};
            
            [~,m_,sp_] = DecodeTensor.special_sess_id_list;<>
            show_filter = ismember(m_, show_mice);
            sp_ = sp_(show_filter);
            m_ = m_(show_filter);
            figure('FileName', 'supplements_pdf/medload/medload_rasters.pdf');
            colorscale = 'log';
            for i = 1:numel(sp_)
                subplot(1,numel(sp_)+1, i);
                t_ = res(sp_(i));
                mean_median_loadings = squeeze(mean(abs(t_.median_loadings)));
                min_d = 30;
                ns = t_.n_sizes;
                im_data = (mean_median_loadings(ns >= min_d,1:min_d));
                %padded_im_data = nan(16 - size(im_data,1), size(im_data,2));
                %imagesc([im_data;padded_im_data], [-1.6 log10(0.3)]);
                surf(1:min_d, ns(ns>=min_d), im_data, 'EdgeColor', 'none');
                view(2);
                
                %set(gca, 'XScale', 'log');
                %set(gca, 'YScale', 'log');
                set(gca, 'ColorScale', colorscale);
                caxis([0.03 0.25]);
                xlim([1 min_d]);
                ylim([min_d+10, 500]);
                xlabel 'Fluctuation mode, i'
                ylabel 'Number of cells'
                title(sprintf('Mouse %s', m_{i}(end-1:end)), 'FontName', 'Helvetica', 'FontSize', 6, 'FontWeight', 'normal', 'Color', 'b');
                
                set(gca, 'FontSize', 6);
                set(gca, 'FontName', 'Helvetica');
                set(gca, 'TickLength', [0.02 0.02]);
                colorbar;
                %set(gca, 'YTick', [1 2 4 8 16 32 64].*min_d);
                %rectangle('Position',...
                %    0.5+[0 size(im_data,1) size(im_data,2) (48 - size(im_data,1))],...
                %    'FaceColor', 'w', 'EdgeColor', 'k', 'LineStyle', 'none');
                %set(gca, 'YTickLabel', 10*cellfun(@str2num, get(gca, 'YTickLabel')));
                box off;
                if i > 1
                    box off
                    xlabel ''
                    ylabel ''
                    set(gca, 'YTick', []);
                end
                %colorbar;
            end
            subplot(1, numel(sp_)+1, numel(sp_)+1);
            t_ = res(sp_(1));
            mean_median_loadings_s = squeeze(mean(abs(t_.median_loadings_s)));
            min_d = 30;
            ns = t_.n_sizes;
            im_data = (mean_median_loadings_s(ns >= min_d,1:min_d));
            surf(1:min_d, ns(ns>=min_d), im_data, 'EdgeColor', 'none');
            view(2);
            %set(gca, 'YScale', 'log');
            set(gca, 'ColorScale', colorscale);
            caxis([0.03 0.25]);
            xlim([1 min_d]);
            ylim([min_d+10, 500]);
            %xlabel 'Fluctuation mode, i'
            %ylabel 'Number of cells'
            title('Shuffled', 'FontName', 'Helvetica', 'FontSize', 6, 'FontWeight', 'normal', 'Color', 'r');
            
            set(gca, 'FontSize', 6);
            set(gca, 'FontName', 'Helvetica');
            set(gca, 'TickLength', [0.02 0.02]);
            set(gca, 'YTick', [1 2 4 8 16 32 64].*min_d);
            box off
            %axis off
            xlabel ''; ylabel '';
            set(gca, 'YTick', []);
            colorbar;
            set(gcf, 'Units', 'inches');
            set(gcf, 'Position', [8.5521    6.2292    8.3125    1.6146]);
            colormap parula;
            Utils.printto;
            
            figure('FileName', 'figure2_pdf/medload/medload_curves_with_mean.pdf');
            MultiSessionVisualizer.plot_single_filtered(n_sizes, series, {'b', 'r'}, sp_);
            set(gca, 'XScale', 'log');
            set(gca, 'YScale', 'log');
            xlabel 'Number of cells'
            ylabel 'mean_i|cos(PC_i, Dm)|, first 5'
            xlim([1 500]);
            ylim([-Inf 1]);
            figure_format([1 1.4]);
            Utils.printto;
            
            MultiSessionVisualizer.plot_series(n_sizes, series, {'b','r'}, mouse_name);
            axs = findall(gcf, 'type', 'axes');
            set(axs, 'YScale', 'log');
            set(axs, 'XScale', 'log');
            xlabel 'Number of cells'
            ylabel 'mean_i|cos(PC_i, Dm)|, first 5'
            multi_figure_format;
            Utils.printto('supplements_pdf/medload', 'multi_medload_curves_with_mean.pdf');
            
            n_c = 50;
            fit_func = @(x,y)fit(log10(x(x>=n_c))',log10(mean(y(:,x>=n_c)))', 'poly1');
            [fr_, gf_] = cellfun(fit_func, n_sizes, series{1}, 'UniformOutput', false);
            [fr_s, gf_s] = cellfun(fit_func, n_sizes, series{2}, 'UniformOutput', false);
            
            gf_ = cell2mat(gf_);
            rsquare = [gf_.rsquare];
            fprintf('Unshuf: range %f-%f, median %f\n', min(rsquare), max(rsquare), median(rsquare));
            
            gf_s = cell2mat(gf_s);
            rsquare_s = [gf_s.rsquare];
            fprintf('Shuf: range %f-%f, median %f\n', min(rsquare_s), max(rsquare_s), median(rsquare_s));
            
            [rate_f, rate_f_conf] = Utils.fit_get(fr_, 'p1');
            [rate_f_s, rate_f_s_conf] = Utils.fit_get(fr_s, 'p1');
            figure('FileName', 'figure2_pdf/medload/inset_with_mean.pdf');
            Utils.bns_groupings(rate_f, rate_f_s, rate_f_conf, rate_f_s_conf, mouse_name, true);
            hold on;
            %line(xlim-0.5, [0 0], 'Color', 'k', 'LineStyle', '-');
            line(xlim, [-0.5 -0.5], 'Color', 'k', 'LineStyle', ':');
            ylabel 'Fit exponent'
            ylim([-Inf Inf]);
            set(gca, 'XTickLabels', {'Real', 'Shuf.'});
            %figure_format([0.8 1]/2, 'fontsize', 4);
            Utils.specific_format('inset');
            Utils.printto;
            
            figure('FileName', 'supplements_pdf/medload/multi_inset_with_mean.pdf');
            Utils.bns_groupings(rate_f, rate_f_s, rate_f_conf, rate_f_s_conf, mouse_name, false, {'Unshuffled', 'Shuffled'});
            hold on;
            line(xlim-0.5, [0 0], 'Color', 'k', 'LineStyle', '-');
            line(xlim, [-0.5 -0.5], 'Color', 'k', 'LineStyle', ':');
            ylabel 'Fit exponent'
            ylim([-Inf Inf]);
            %set(gca, 'XTickLabels', {'Unsh.', 'Sh.'});
            %figure_format([0.8 1]/2, 'fontsize', 4);
            Utils.specific_format('MBNS');
            Utils.printto;
            %}
            %PABLO POINT 4
            if true
                fprintf('The following are the values (there is no variation for max cells)\n for max_i cos(PC_i, Dm), at the maximal # of cells\n');
                max_overlap = cellfun(@(x)mean(x(:,end)), series{1}); %vector of length 107
                disp(max_overlap);
                
                
                fit_savefile = 'decoding_curves_fits.mat';
                recompute = false;
                if recompute || ~exist(fit_savefile, 'file')
                    PanelGenerator.decoding_curves('remake', true, 'recompute', recompute);
                end
                load(fit_savefile);
                
                fprintf('The following is Pablo''s formula: (I0N(shuf)-I0N(real))\n');
                pablos_formula = (I0_fit_s.*N_fit_s - I0_fit.*N_fit) ./ (I0_fit_s.*N_fit_s + I0_fit.*N_fit);
                disp(pablos_formula);
                %{
                figure;
                scatter(max_overlap, pablos_formula);
                xlabel 'mean_i|cos(PC_i, Dm)|, first 5'
                ylabel 'I_0N_{sh}-I_0N_{re}'
                refline
                [res_, GOF_] = fit(max_overlap(:), pablos_formula(:), 'poly1');
                text(0.25, -0.6, sprintf('\\it{R}^2 = %.2e', GOF_.rsquare));
                %}
                [~, mouse_names] = DecodeTensor.filt_sess_id_list;
                figure;
                PanelGenerator.plot_regress(max_overlap, pablos_formula,...
                    0*max_overlap, 0*pablos_formula, mouse_names, 'g', 'text_coords', [0.15 -0.5]);
                xlabel 'mean_i|cos(PC_i, \Delta\mu)|, first 5'
                ylabel '(I_0N_{sh}-I_0N_{re}) / (I_0N_{sh}+I_0N_{re})'
                figure_format('factor', 1.6);
                Utils.printto('supplements_pdf', 'Pablo_point_4_mean_RATIO.pdf');
                
                
                figure;
                PanelGenerator.plot_regress_averaged(max_overlap, pablos_formula,...
                    0*max_overlap, 0*pablos_formula, mouse_names, 'g', 'text_coords', [0.15 -0.5]);
                xlabel 'mean_i|cos(PC_i, \Delta\mu)|, first 5'
                ylabel '(I_0N_{sh}-I_0N_{re}) / (I_0N_{sh}+I_0N_{re})'
                figure_format('factor', 1.6);
                Utils.printto('supplements_pdf', 'Pablo_point_4_mean_grouped_RATIO.pdf');
            end            
        end
        
        function m = mean_func(x)
            x(x==0) = nan;   
            m = nanmean(x,3); 
        end
        
        function m = mean_func_trunc5(x)
            x(x==0) = nan;   
            m = nanmean(x(:,:,1:5),3); 
        end
        %function adjacent_decoders
        %    load('adjacent_agg_190725-094031_0.mat');
        %    keyboard
        %    n_sizes = DecodeTensor.get_n_neurons_filt([res.filt_index]);
        %    n_sizes = arrayfun(@(x)[2, 10:10:x, x], n_sizes, 'UniformOutput', false);
        %    m2_slope = Utils.cf_(@Utils.fitaline, n_sizes, {res.m2});
        %    m2_s_slope = Utils.cf_(@Utils.fitaline, n_sizes, {res.m2_s});
        %
        %    n_cutoff = 100;
        %    nv_slope = Utils.cf_(@(n,m)Utils.fitaline(n,m,n_cutoff), n_sizes, {res.nv});
        %    nv_s_slope = Utils.cf_(@(n,m)Utils.fitaline(n,m,n_cutoff), n_sizes, {res.nv_s});
        %
        %    asymp_dp2_w = Utils.cf_(@(x,y)x./y, m2_slope, nv_slope);
        %    asymp_dp2_w_s = Utils.cf_(@(x,y)x./y, m2_s_slope, nv_s_slope);
        %end
        
        function adjacent_decoders(recompute)
            if ~exist('recompute', 'var')
                recompute = false;
            end
            load('adjacent_agg_190725-194911_0.mat');
            %keyboard
            n_cutoff = 100;
            medify = @(z) Utils.cf_(@(y)cellfun(@(x)median(x(:)), y), z);
            %medify = @(z) Utils.cf_(@(y)cellfun(@(x)mean(x(:)), y), z);
            full_line = @Utils.fitaline;
            asymp_line = @(n,m) Utils.fitaline(n,m,n_cutoff);
            asymp_line_intercept = @(n,m) Utils.fitaline(n,m,n_cutoff,true);
            
            [m2_slopes, m2_slope_confs] = cellfun(full_line, {res.n_sizes}, medify({res.m2}), 'UniformOutput', false);
            [m2_s_slopes, m2_s_slope_confs] = cellfun(full_line, {res.n_sizes}, medify({res.m2_s}), 'UniformOutput', false);
            [m2_d_slopes, m2_d_slope_confs] = cellfun(full_line, {res.n_sizes}, medify({res.m2_d}), 'UniformOutput', false);
            
            [nv_slopes, nv_slope_confs] = cellfun(asymp_line, {res.n_sizes}, medify({res.nv}), 'UniformOutput', false);
            [nv_s_slopes, nv_s_slope_confs] = cellfun(asymp_line, {res.n_sizes}, medify({res.nv_s}), 'UniformOutput', false);
            [nv_d_slopes, nv_d_slope_confs] = cellfun(asymp_line, {res.n_sizes}, medify({res.nv_d}), 'UniformOutput', false);
            [nv_intercepts, nv_intercept_confs] = cellfun(asymp_line_intercept, {res.n_sizes}, medify({res.nv}), 'UniformOutput', false);
            [nv_s_intercepts, nv_s_intercept_confs] = cellfun(asymp_line_intercept, {res.n_sizes}, medify({res.nv_s}), 'UniformOutput', false);
            [nv_d_intercepts, nv_d_intercept_confs] = cellfun(asymp_line_intercept, {res.n_sizes}, medify({res.nv_d}), 'UniformOutput', false);
            
            dp2_w = cellfun(@(x,y)x/y, m2_slopes, nv_slopes);
            dp2_w_conf = cellfun(@(x, xc, y, yc) abs(x/y)*sqrt((xc/x)^2 + (yc/y)^2), m2_slopes, m2_slope_confs, nv_slopes, nv_slope_confs);
            dp2_ws = cellfun(@(x,y)x/y, m2_s_slopes, nv_s_slopes);
            dp2_ws_conf = cellfun(@(x, xc, y, yc) abs(x/y)*sqrt((xc/x)^2 + (yc/y)^2), m2_s_slopes, m2_s_slope_confs, nv_s_slopes, nv_s_slope_confs);
            dp2_wd = cellfun(@(x,y)x/y, m2_d_slopes, nv_d_slopes);
            dp2_wd_conf = cellfun(@(x, xc, y, yc) abs(x/y)*sqrt((xc/x)^2 + (yc/y)^2), m2_d_slopes, m2_d_slope_confs, nv_d_slopes, nv_d_slope_confs);
            
            dp2_w_slope = cellfun(@(x,y)x/y, m2_slopes, nv_intercepts);
            dp2_w_slope_conf = cellfun(@(x, xc, y, yc) abs(x/y)*sqrt((xc/x)^2 + (yc/y)^2), m2_slopes, m2_slope_confs, nv_intercepts, nv_intercept_confs);
            dp2_ws_slope = cellfun(@(x,y)x/y, m2_s_slopes, nv_s_intercepts);
            dp2_ws_slope_conf = cellfun(@(x, xc, y, yc) abs(x/y)*sqrt((xc/x)^2 + (yc/y)^2), m2_s_slopes, m2_s_slope_confs, nv_s_intercepts, nv_s_intercept_confs);
            dp2_wd_slope = cellfun(@(x,y)x/y, m2_d_slopes, nv_d_intercepts);
            dp2_wd_slope_conf = cellfun(@(x, xc, y, yc) abs(x/y)*sqrt((xc/x)^2 + (yc/y)^2), m2_d_slopes, m2_d_slope_confs, nv_d_intercepts, nv_d_intercept_confs);
            
            
            fit_savefile = 'decoding_curves_fits.mat';
            if recompute || ~exist(fit_savefile, 'file')
                PanelGenerator.decoding_curves('remake', true, 'recompute', recompute);
            end
            load(fit_savefile);
            
            good_fit_filter = (I0_conf < 0.5*I0_fit) &...
                (I0_conf_s < 0.5*I0_fit_s) &...
                (N_conf < 0.5*N_fit);% & (N_fit < 500) & (cellfun(@max, {res.n_sizes}) >= 300); %200
            g_ = good_fit_filter;
            %g_ = true(size(g_));
            
            %figure;
            InfoLimit = N_fit.*I0_fit;
            InfoLimit_conf = abs(InfoLimit).*sqrt((N_conf./N_fit).^2 + (I0_conf./I0_fit).^2);
            
            x_cell = Utils.cf_(@(x)x(g_), {dp2_w, dp2_ws, dp2_wd, dp2_w_slope, dp2_ws_slope, dp2_wd_slope});
            x_conf_cell = Utils.cf_(@(x)x(g_), {dp2_w_conf, dp2_ws_conf, dp2_wd_conf, dp2_w_slope_conf, dp2_ws_slope_conf, dp2_wd_slope_conf});
            x_lab_cell = {'dp2_w', 'dp2_{ws}', 'dp2_{wd}', 'dp2_w slope', 'dp2_{ws} slope', 'dp2_{wd} slope'};
            
            y_cell = Utils.cf_(@(x)x(g_), {I0_fit, N_fit, InfoLimit});
            y_conf_cell = Utils.cf_(@(x)x(g_), {I0_conf, N_conf, InfoLimit_conf});
            y_lab_cell = {'I_0', 'N', 'I_0N'};
            
            %PanelGenerator.aux_many_regress(x_cell(1:2:end), y_cell(1:2:end), x_conf_cell(1:2:end), y_conf_cell(1:2:end), mouse_names(g_), x_lab_cell(1:2:end), y_lab_cell(1:2:end));
            PanelGenerator.aux_many_regress(x_cell, y_cell, x_conf_cell, y_conf_cell, mouse_names(g_), x_lab_cell, y_lab_cell);
            %keyboard;
            
            figure; PanelGenerator.plot_regress(I0_fit(g_), InfoLimit(g_), I0_conf(g_), InfoLimit_conf(g_), mouse_names(g_), 'k');
            xlabel I_0; ylabel I_0N
            
            figure; PanelGenerator.plot_regress(dp2_wd(g_), dp2_w(g_), dp2_wd_conf(g_), dp2_w_conf(g_), mouse_names(g_), 'k');
            xlabel 'dp2_{wd}'; ylabel 'dp2_w';
            fprintf('Using n = %d samples\n', sum(g_));
            
            

            
            
            
            figure('FileName', 'figure2_pdf/signal_and_noise/noise_curves_adjacent.pdf');
            
            show_mice = {'Mouse2022', 'Mouse2024', 'Mouse2028'};
            
            %[~,m_,sp_] = DecodeTensor.special_sess_id_list;<>
            %show_filter = ismember(m_, show_mice);
            %sp_ = sp_(show_filter);
            sp_ = SessManager.special_sessions_usable_index(show_mice);
            %[md2_slopes, ~] = cellfun(full_line, {res.n_sizes}, medify({res.md2}), 'UniformOutput', false);
            n_f_ = @(x) Utils.cf_(@(y,z)y./z, medify(x), m2_d_slopes);
            MultiSessionVisualizer.plot_single_filtered({res.n_sizes}, {n_f_({res.m2_d}), n_f_({res.nv_d}), n_f_({res.nv_s})}, {'k', 'b', 'r'}, sp_);
            xlabel 'Number of cells'
            
            ylabel(sprintf('s^2/K (num. of cells)'));
            
            text(20, 370/2, '(Dm)^2', 'Color', 'k', 'HorizontalAlignment', 'left');
            text(20, 470/2, 's^2', 'Color', 'b', 'HorizontalAlignment', 'left');
            text(20, 570/2, 's^2', 'Color', 'r', 'HorizontalAlignment', 'left');
            %figure_format('factor', 1.3);
            figure_format([1 1.4], 'factor', 1);
            Utils.printto;
            

            
            figure('FileName', 'figure2_pdf/signal_and_noise/noise_rate_of_change_adjacent.pdf');
            kick_out = ismember(mouse_names, {'Mouse2010', 'Mouse2011', 'Mouse2012', 'Mouse2021'});
            kick_out = kick_out & false;
            Utils.bns_groupings(cell2mat(nv_d_slopes(~kick_out))./cell2mat(m2_d_slopes(~kick_out)), cell2mat(nv_s_slopes(~kick_out))./cell2mat(m2_d_slopes(~kick_out)), cell2mat(nv_d_slope_confs(~kick_out))./cell2mat(m2_d_slopes(~kick_out)), cell2mat(nv_s_slope_confs(~kick_out))./cell2mat(m2_d_slopes(~kick_out)), mouse_names(~kick_out), true);
            ylim([-Inf Inf]);
            ylabel(sprintf('s^2/K slope'));
            figure_format([0.666 1.4]);
            Utils.printto;
            
       
            figure('FileName', 'supplements_pdf/signal_and_noise/multi_noise_rate_of_change_adjacent.pdf');
            Utils.bns_groupings(cell2mat(nv_d_slopes)./cell2mat(m2_d_slopes), cell2mat(nv_s_slopes)./cell2mat(m2_d_slopes), cell2mat(nv_d_slope_confs)./cell2mat(m2_d_slopes), cell2mat(nv_s_slope_confs)./cell2mat(m2_d_slopes), mouse_names, false);
            ylim([-Inf Inf]);
            ylabel(sprintf('s^2 rate of change'));
            Utils.specific_format('MBNS');
            Utils.printto;
            
            

            PanelGenerator.aux_regressions('limit_vs_dp2_wd.pdf', dp2_wd(g_), InfoLimit(g_),...
                dp2_wd_conf(g_), InfoLimit_conf(g_), mouse_names(g_), 'b',...
                [1 0.08], '(Dm)^2 slope / s^2 slope', 'I_0N (cm^{-2})',...
                false, 'figure2_pdf/adjacent', 'supplements_pdf/adjacent');
            
            PanelGenerator.aux_regressions('limit_vs_dp2_w.pdf', dp2_w(g_), InfoLimit(g_),...
                dp2_w_conf(g_), InfoLimit_conf(g_), mouse_names(g_), 'b',...
                [1 0.08], '(Dm)^2 slope / s^2 slope (on w)', 'I_0N (cm^{-2})',...
                false, 'figure2_pdf/adjacent', 'supplements_pdf/adjacent');
            
            PanelGenerator.aux_regressions('I0_vs_dp2_ws_slope.pdf', dp2_ws_slope(g_), I0_fit(g_),...
                dp2_ws_slope_conf(g_), I0_conf(g_), mouse_names(g_), 'r',...
                [0.01 6e-4], '(Dm_s)^2 slope / s_s^2 intercept',...
                sprintf('I_0 fit value\n(cm^{-2}neuron^{-1})'), 'y', 'figure2_pdf/adjacent',...
                'supplements_pdf/adjacent');
            
            PanelGenerator.aux_regressions('I0_vs_limit.pdf', I0_fit(g_), InfoLimit(g_),...
                I0_conf(g_), InfoLimit_conf(g_), mouse_names(g_), 'k',...
                [6e-4 0.08], sprintf('I_0 fit value\n(cm^{-2}neuron^{-1})'),...
                'I_0N (cm^{-2})', 'x', 'figure2_pdf/adjacent', 'supplements_pdf/adjacent');
            
            PanelGenerator.aux_regressions('dp2_ws_slope_vs_dp2_wd.pdf', dp2_ws_slope(g_),...
                dp2_wd(g_), dp2_ws_slope_conf(g_), dp2_wd_conf(g_), mouse_names(g_), 'k',...
                [0.01 1], '(Dm_s)^2 slope / s_s^2 intercept',...
                '(Dm)^2 slope / s^2 slope', false, 'figure2_pdf/adjacent', 'supplements_pdf/adjacent');
            
            PanelGenerator.aux_regressions('dp2_wd_vs_dp2_w.pdf', dp2_wd(g_),...
                dp2_w(g_), dp2_wd_conf(g_), dp2_w_conf(g_), mouse_names(g_), 'k',...
                [0.01 1], '(Dm)^2 slope / s^2 slope',...
                '(Dm)^2 slope / s^2 slope (on w)', false, 'figure2_pdf/adjacent', 'supplements_pdf/adjacent');
            
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
            
            if filt_num
                figure('FileName', 'figure2_pdf/signal_and_noise_filt_num/noise_rate_of_change.pdf');
            else
                figure('FileName', 'figure2_pdf/signal_and_noise/noise_rate_of_change.pdf');
            end
            kick_out = ismember(mouse_list, {'Mouse2010', 'Mouse2011', 'Mouse2012', 'Mouse2021'});
            Utils.bns_groupings(sm_slope(~kick_out), sms_slope(~kick_out), sm_conf(~kick_out), sms_conf(~kick_out), mouse_list(~kick_out), true);
            ylim([-Inf Inf]);
            ylabel(sprintf('s^2 along Dm\nrate of change'));
            figure_format('factor', 1.6);
            Utils.printto;
            
            if filt_num
                figure('FileName', 'supplements_pdf/signal_and_noise_filt_num/multi_noise_rate_of_change.pdf');
            else
                figure('FileName', 'supplements_pdf/signal_and_noise/multi_noise_rate_of_change.pdf');
            end
            Utils.bns_groupings(sm_slope, sms_slope, sm_conf, sms_conf, mouse_list, false);
            ylim([-Inf Inf]);
            ylabel(sprintf('s^2 along Dm\nrate of change'));
            Utils.specific_format('MBNS');
            Utils.printto;
            %Utils.create_svg(gcf, 'figure2_svg', 'grouped_noise_rate_of_change');
            
            if filt_num
                figure('FileName', 'figure2_pdf/signal_and_noise_filt_num/noise_curves.pdf');
            else
                figure('FileName', 'figure2_pdf/signal_and_noise/noise_curves.pdf');
            end
            show_mice = {'Mouse2022', 'Mouse2024', 'Mouse2028'};
            
            %[~,m_,sp_] = DecodeTensor.special_sess_id_list;<>
            %show_filter = ismember(m_, show_mice);
            %sp_ = sp_(show_filter);
            sp_ = SessManager.special_sessions_usable_index(show_mice);
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
            
            good_fit_filter = (I0_upper < 0.5*I0_fit_value) &...
                (I0_upper_s < 0.5*I0_fit_value_s) &...
                (N_upper < 0.5*N_fit_value);% & (N_fit < 500) & (cellfun(@max, {res.n_sizes}) >= 300); %200
            g_ = good_fit_filter;
            disp(find(~g_));
            if filt_num
                fname = 'imse_limit_regression_filtnum.pdf';
            else
                fname = 'imse_limit_regression.pdf';
            end
            %scatter(I0_fit_value.*N_fit_value, 1./sm_slope, 4, 'b');
            limit_uncertainty = sqrt((I0_fit_value.*(N_upper)).^2 + (N_fit_value.*(I0_upper)).^2);
            inv_sm_slope_uncertainty = sm_conf ./ sm_slope.^2;
            hold on;
            
            %errorbar(I0_fit_value(g_).*N_fit_value(g_), 1./sm_slope(g_), inv_sm_slope_uncertainty(g_), inv_sm_slope_uncertainty(g_),...
            %    limit_uncertainty(g_), limit_uncertainty(g_), 'LineStyle', 'none', 'Color', 'k', 'CapSize', 1);
            %scatter(I0_fit_value(g_).*N_fit_value(g_), 1./sm_slope(g_), dotsize, DecodeTensor.mcolor(mouse_names(g_), false), 'filled');
            %[fitresult, adjr2] = Utils.regress_line(I0_fit_value(g_).*N_fit_value(g_), 1./sm_slope(g_));
            %h_ = plot(fitresult); legend off
            %h_.Color = 'b';
            %xlim([-Inf 0.15]);
            %text(0.1, 1, sprintf('adj. R^2 = %.2f', adjr2));
            %xlabel 'IMSE limit I_0N';
            %ylabel(sprintf('Inverse s^2\nrate of change'));
            %fprintf('imse limit regression, N = %d\n', numel(I0_fit_value(g_)));
            %figure_format('factor', 1.6);
            %%Utils.create_svg(gcf, 'figure2_svg', 'imse_limit_regression');
            %Utils.printto;
            PanelGenerator.aux_regressions(fname, 1./sm_slope(g_), I0_fit_value(g_).*N_fit_value(g_),...
                inv_sm_slope_uncertainty(g_), limit_uncertainty(g_), mouse_names(g_), 'b',...
                [4 0.08], sprintf('Inverse s^2/K\nrate of change'), 'I_0N (cm^{-2})',...
                false, 'figure2_pdf/signal_and_noise', 'supplements_pdf/signal_and_noise');
            
            
            if filt_num
                fname = 'I0_value_regression_filtnum.pdf';
            else
                fname = 'I0_value_regression.pdf';
            end
            %scatter(I0_fit_value_s, 1./sms_intercept, 'r');
            %hold on;
            inv_intercept_errb = sms_intercept_conf./sms_intercept.^2;
            PanelGenerator.aux_regressions(fname, 1./sms_intercept(g_), I0_fit_value(g_), inv_intercept_errb(g_), I0_upper(g_), mouse_names(g_), 'r',...
                [0 6e-4], sprintf('Inverse s^2/K intercept'), sprintf('I_0 fit value\n(cm^{-2}neuron^{-1})'), 'y', 'figure2_pdf/signal_and_noise', 'supplements_pdf/signal_and_noise');
            
            %errorbar(I0_fit_value_s(g_), 1./sms_intercept(g_), inv_intercept_errb(g_), inv_intercept_errb(g_), I0_upper_s(g_), I0_upper_s(g_), 'LineStyle', 'none', 'Color', 'k', 'CapSize', 1);
            %scatter(I0_fit_value_s(g_), 1./sms_intercept(g_), dotsize, DecodeTensor.mcolor(mouse_names(g_), false), 'filled');
            %[fitresult, adjr2] = Utils.regress_line(I0_fit_value_s(g_), 1./sms_intercept(g_));
            %plot(fitresult); legend off
            %text(7e-4, 0.055, sprintf('adj. R^2 = %.2f', adjr2));
            %xlabel 'I_0 fit value'
            %ylabel 'Asymptotic 1/s^2'
            %set(gca, 'XTickLabel', arrayfun(@Utils.my_fmt, get(gca, 'XTick') ,'UniformOutput', false));
            %fprintf('I0 value regression, N = %d\n', numel(I0_fit_value_s(g_)));
            %figure_format('factor', 1.6);
            %%Utils.create_svg(gcf, 'figure2_svg', 'I0_value_regression');
            %Utils.printto;
        end
    end
end