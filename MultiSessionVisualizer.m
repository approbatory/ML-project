%%TODO: this class will allow visualization of all sessions,
%by mouse, for whichever quantities desired, of the raw curve
%rather than the parameter fit

classdef MultiSessionVisualizer
    methods(Static)
        function SnN(normed, filt_num)
            load('signal_and_noise_final.mat');
            [~, mouse_list] = DecodeTensor.filt_sess_id_list;
            
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
                ylabel(sprintf('Distance^2\n(on \DeltaF/F values)'));
            end
        end
    
        function plot_series(n_sizes, series_cell, color_cell, mouse_list)
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
                        %if n_{k}(end)<200 %%FILTER FOR <200
                        %    continue;
                        %end
                        Utils.neuseries(n_{k}, s_{k}, c_);
                        hold on;
                    end %session in mouse
                end %quantity shown
                xlim([0 500]);
                ylim([0 Inf]);
                title(mouse_names{m_i});
            end %mice id
        end %func
    
    end
end
