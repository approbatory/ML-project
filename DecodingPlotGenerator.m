classdef DecodingPlotGenerator < handle
    properties
        dbfile;
        figure_collection;
    end
    methods
        function o = DecodingPlotGenerator(dbfile)
            o.dbfile = dbfile;
            o.figure_collection = [];
        end
        
        function pairwise_confidence_plots(o, mouse_name)
            conn = sqlite(o.dbfile);
            cleaner = onCleanup(@()conn.close);
            
            command_template = @(isbwA, isbwB, offset, setting, mouse)...
                sprintf(['select avg(CorrectConfidence) from pairwise'...
                ' where BinA between %d*20+1 and %d*20+20 and'...
                ' BinB between %d*20+1 and %d*20+20 and'...
                ' BinB = BinA + %d*20 + %d and'...
                ' Setting = ''%s'' and Mouse = ''%s'' group by BinA, BinB;'],...
                isbwA, isbwA, isbwB, isbwB, isbwB && ~isbwA, offset, setting, mouse);
            
            same_conf = zeros(2,20); same_conf_err = zeros(2,20);
            diff_conf = zeros(2,20); diff_conf_err = zeros(2,20);
            settings = {'unshuffled', 'shuffled'};
            for setting_code = 1:2
                setting = settings{setting_code};
                for offset = 0:19
                    ret_right = cell2mat(conn.fetch(command_template(0, 0, offset, setting, mouse_name)));
                    ret_left = cell2mat(conn.fetch(command_template(1, 1, offset, setting, mouse_name)));
                    ret = [ret_right ; ret_left];
                    same_conf(setting_code, offset+1) = mean(ret);
                    same_conf_err(setting_code, offset+1) = std(ret) ./ sqrt(length(ret));
                    
                    ret_diff = cell2mat(conn.fetch(command_template(0, 1, offset, setting, mouse_name)));
                    diff_conf(setting_code, offset+1) = mean(ret_diff);
                    diff_conf_err(setting_code, offset+1) = std(ret_diff) ./ sqrt(length(ret_diff));
                end
            end
            
            o.figure_collection = [o.figure_collection, figure('FileName', 'pairwise_confidence_by_distance.png')];
            hold on;
            errorbar(0:19, same_conf(1,:), same_conf_err(1,:), '-b', 'DisplayName', 'Same direction');
            errorbar(0:19, same_conf(2,:), same_conf_err(2,:), '-r', 'DisplayName', 'Same direction, shuffled');
            
            errorbar(0:19, diff_conf(1,:), diff_conf_err(1,:), ':b', 'DisplayName', 'Opposite directions');
            errorbar(0:19, diff_conf(2,:), diff_conf_err(2,:), ':r', 'DisplayName', 'Opposite directions, shuffled');
            ylim([0.5 1]);
            legend; legend boxoff;
            
            xlabel 'Bin distance'
            ylabel 'Correct bin posterior'
            
            title(sprintf('Confidence level of correct bin in pairwise SVM decoding,\n from %s', mouse_name));
        end
        
        function series_NumNeurons(o, varargin) %controls for datasize by default
            p = inputParser;
            p.addOptional('control_DataSize', 3, @isscalar); %options are 3 or false (for no control)
            p.parse(varargin{:});
            namer = @(n) sprintf('%s_NumNeurons_%d.png', n, p.Results.control_DataSize);
            
            conn = sqlite(o.dbfile);
            cleaner = onCleanup(@()conn.close);
            
            if p.Results.control_DataSize == false
                dataSize_param = 'max';
                Mouse_values = conn.fetch('select distinct Mouse from decoding order by Mouse;');
            else
                command = 'select max(DataSize) from decoding group by Mouse';
                data_sizes = cell2mat(conn.fetch(command));
                dataSize_param = min(maxk(data_sizes, p.Results.control_DataSize));
                command = sprintf('select distinct Mouse from decoding where DataSize >= %d order by Mouse;', dataSize_param);
                Mouse_values = conn.fetch(command);
            end
            
            %plotting neuron series mean error shuf/unshuf
            o.figure_collection = [o.figure_collection, figure('FileName', namer('mean_errs'))];
            DecodingPlotGenerator.plot_mice(conn, Mouse_values, 'shuffled', 'NumNeurons', 'MeanErrors', dataSize_param);
            DecodingPlotGenerator.plot_mice(conn, Mouse_values, 'unshuffled', 'NumNeurons', 'MeanErrors', dataSize_param);
            %decorations
            set(gca, 'YScale', 'log');
            xlabel 'Number of cells'
            ylabel 'Mean error (cm)'
            xlim([0 500]);
            text(200, 7, 'Unshuffled', 'Color', 'blue');
            text(300, 2.2, 'Shuffled', 'Color', 'red');
            %end decorations
            
            %plotting neuron series mean error diag/unshuf
            o.figure_collection = [o.figure_collection, figure('FileName', namer('diag_mean_errs'))];
            DecodingPlotGenerator.plot_mice(conn, Mouse_values, 'diagonal', 'NumNeurons', 'MeanErrors', dataSize_param);
            DecodingPlotGenerator.plot_mice(conn, Mouse_values, 'unshuffled', 'NumNeurons', 'MeanErrors', dataSize_param);
            %decorations
            set(gca, 'YScale', 'log');
            xlabel 'Number of cells'
            ylabel 'Mean error (cm)'
            xlim([0 500]); ylim([1 100]);
            text(200, 8, 'Diagonal', 'Color', 'magenta');
            text(300, 2.2, 'Full', 'Color', 'blue');
            %end decorations
            
            %plotting neuron series IMSE shuf/unshuf
            o.figure_collection = [o.figure_collection, figure('FileName', namer('IMSE'))];
            DecodingPlotGenerator.plot_mice(conn, Mouse_values, 'shuffled', 'NumNeurons', 'IMSE', dataSize_param);
            DecodingPlotGenerator.plot_mice(conn, Mouse_values, 'unshuffled', 'NumNeurons', 'IMSE', dataSize_param);
            %decorations
            xlabel 'Number of cells'
            ylabel '1/MSE (cm^{-2})'
            xlim([0 500]);
            text(300, 0.2, 'Unshuffled', 'Color', 'blue');
            text(400, 0.6, 'Shuffled', 'Color', 'red');
            %end decorations
        end
        
        
        function series_DataSize(o, varargin)
            p = inputParser;
            p.addOptional('control_NumNeurons', false, @isscalar);
            p.parse(varargin{:});
            namer = @(n) sprintf('%s_DataSize_%d.png', n, p.Results.control_NumNeurons);
            
            conn = sqlite(o.dbfile);
            cleaner = onCleanup(@()conn.close);
            
            if p.Results.control_NumNeurons == false
                numNeurons_param = 'max';
                Mouse_values = conn.fetch('select distinct Mouse from decoding where DataSize = 2 order by Mouse;');
            else
                error('not supported');
                %command = sprintf('select min(DS) from (select max(DataSize) as DS from decoding group by Mouse order by max(DataSize) desc limit %d);', p.Results.control_DataSize);
                command = 'select max(NumNeurons) from decoding group by Mouse';
                num_neurons = cell2mat(conn.fetch(command));
                numNeurons_param = min(maxk(num_neurons, p.Results.control_NumNeurons));
                command = sprintf('select distinct Mouse from decoding where NumNeurons >= %d;', numNeurons_param);
                Mouse_values = conn.fetch(command);
            end
            
            %plotting neuron series mean error shuf/unshuf
            o.figure_collection = [o.figure_collection, figure('FileName', namer('mean_errs'))];
            DecodingPlotGenerator.plot_mice(conn, Mouse_values, 'shuffled', 'DataSize', 'MeanErrors', numNeurons_param);
            DecodingPlotGenerator.plot_mice(conn, Mouse_values, 'unshuffled', 'DataSize', 'MeanErrors', numNeurons_param);
            %decorations
            set(gca, 'YScale', 'log');
            xlabel 'Number of trials'
            ylabel 'Mean error (cm)'
            %xlim([0 500]);
            text(100, 3.5, 'Unshuffled', 'Color', 'blue');
            text(100, 2, 'Shuffled', 'Color', 'red');
            %end decorations
            
            %plotting neuron series mean error diag/unshuf
            o.figure_collection = [o.figure_collection, figure('FileName', namer('diag_mean_errs'))];
            DecodingPlotGenerator.plot_mice(conn, Mouse_values, 'diagonal', 'DataSize', 'MeanErrors', numNeurons_param);
            DecodingPlotGenerator.plot_mice(conn, Mouse_values, 'unshuffled', 'DataSize', 'MeanErrors', numNeurons_param);
            %decorations
            set(gca, 'YScale', 'log');
            xlabel 'Number of trials'
            ylabel 'Mean error (cm)'
            %xlim([0 500]); ylim([1 100]);
            text(100, 5, 'Diagonal', 'Color', 'magenta');
            text(125, 3, 'Full', 'Color', 'blue');
            %end decorations
            
            %plotting neuron series IMSE shuf/unshuf
            o.figure_collection = [o.figure_collection, figure('FileName', namer('IMSE'))];
            DecodingPlotGenerator.plot_mice(conn, Mouse_values, 'shuffled', 'DataSize', 'IMSE', numNeurons_param);
            DecodingPlotGenerator.plot_mice(conn, Mouse_values, 'unshuffled', 'DataSize', 'IMSE', numNeurons_param);
            %decorations
            xlabel 'Number of trials'
            ylabel '1/MSE (cm^{-2})'
            %xlim([0 500]);
            text(100, 0.3, 'Unshuffled', 'Color', 'blue');
            text(100, 0.6, 'Shuffled', 'Color', 'red');
            %end decorations
        end
        
        function save_figs(o, save_dir)
            if ~exist(save_dir, 'dir')
                mkdir(save_dir);
            end
            for i = 1:numel(o.figure_collection)
                f = o.figure_collection(i);
                figure(f);
                print('-dpng', '-r300', fullfile(save_dir, f.FileName));
            end
        end
        
    end
    
    methods(Static)
        function plot_mice(conn, Mouse_values, setting, series, error_type, uniform_param)
            handles = [];
            for i = 1:numel(Mouse_values)
                mouse = Mouse_values{i};
                [n,m,e] = DecodingPlotGenerator.get_errors(series, conn, mouse, setting, error_type, uniform_param);
                h = DecodingPlotGenerator.errors_plotter(n,m,e,setting, 'index', i, 'DisplayName', mouse);
                handles = [handles h];
            end
            legend(handles);
            legend Location east
            legend boxoff
        end
        function success = errors_plotter(n,m,e,setting, varargin)
            p = inputParser;
            p.addOptional('confidence', 0.95, @isscalar);
            p.addOptional('index', 0, @isscalar);
            p.addOptional('DisplayName', '', @ischar);
            p.parse(varargin{:});

            switch setting
                case 'shuffled'
                    color = 'r';
                case 'unshuffled'
                    color = 'b';
                case 'diagonal'
                    color = 'm';
            end
            assert(length(n) == length(e), 'equal size');
            if length(n) <= 1
                success = [];
                return;
            end
            code = 'ox*sd>ph';
            if p.Results.index == 0
                lp = color;
            else
                lp = [color '-' code(mod(p.Results.index-1, length(code))+1)];
            end
            l_ = shadedErrorBar(n,m,e.*norminv((1+p.Results.confidence)/2),...
                'lineprops', lp);
            l_.mainLine.DisplayName = p.Results.DisplayName;
            success = l_.mainLine;
            %legend(success);
        end
        function [series_values, error_m, error_e] =...
                get_errors(series_type, conn, mouse, setting, error_type, uniform_param)
            if strcmp(error_type, 'IMSE')
                error_type = 'MSE';
                special = 'inverse';
            end
            switch series_type
                case 'NumNeurons'
                    res = conn.fetch(DecodeTensor.build_command(mouse, setting, error_type, [], uniform_param));
                    NumNeurons = cell2mat(res(:,1));
                    DataSize = cell2mat(res(:,2));
                    Error = cell2mat(res(:,3));
                    Series = NumNeurons;
                    Uniform = DataSize;
                    Uniform_name = 'DataSize';
                case 'DataSize'
                    res = conn.fetch(DecodeTensor.build_command(mouse, setting, error_type, uniform_param, []));
                    NumNeurons = cell2mat(res(:,1));
                    DataSize = cell2mat(res(:,2));
                    Error = cell2mat(res(:,3));
                    Series = DataSize;
                    Uniform = NumNeurons;
                    Uniform_name = 'NumNeurons';
                otherwise
                    error('Only NumNeurons or DataSize supported as series_type');
            end
            %switch series_type
            %    case 'NumNeurons'
            %        Series = NumNeurons;
            %        Uniform = DataSize;
            %        Uniform_name = 'DataSize';
            %    case 'DataSize'
            %        Series = DataSize;
            %        Uniform = NumNeurons;
            %        Uniform_name = 'NumNeurons';
            %    otherwise
            %            error('Only NumNeurons or DataSize supported as series_type');
            %end
            assert(length(Error) > 1, 'series unavailable');
            assert(all(Uniform(1) == Uniform), '%s must be uniform', Uniform_name);
            
            series_values = unique(Series);
            error_m = zeros(size(series_values));
            error_e = zeros(size(series_values));
            for i = 1:numel(series_values)
                error_samples = Error(Series == series_values(i));
                if exist('special', 'var') && strcmp(special, 'inverse')
                    error_samples = 1./error_samples;
                end
                error_m(i) = mean(error_samples);
                error_e(i) = std(error_samples) ./ sqrt(length(error_samples));
            end
        end
    end

end