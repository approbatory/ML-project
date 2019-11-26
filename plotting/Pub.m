classdef Pub < handle
    properties
        Q;
        h;
        xlabs;
        ylabs;
        titles;
        subplots;
    end
    methods
        function o = Pub(varargin)
            p = inputParser;
            p.addRequired('width', @isscalar);
            p.addRequired('height', @isscalar);
            p.addParameter('units', 'centimeters', @ischar);
            p.addParameter('hmarginl', 0.08, @isscalar);
            p.addParameter('hmarginr', 0.02, @isscalar);
            p.addParameter('vmargint', 0.09, @isscalar);
            p.addParameter('vmarginb', 0.1, @isscalar);
            p.addParameter('hspacing', 0.11, @isscalar);
            p.addParameter('vspacing', 0.20, @isscalar);
            p.addParameter('layout', []);
            p.addParameter('rows', 1, @isscalar);
            p.addParameter('columns', 1, @isscalar);
            p.parse(varargin{:});
            r = p.Results;
            o.Q = ComputeSubplotPositions(r.rows, r.columns, r.layout,...
                r.hmarginl, r.hmarginr, r.hspacing,...
                r.vmargint, r.vmarginb, r.vspacing);
            o.h = figure;
            ResizeFigure(o.h, r.width, r.height, r.units);
            set(o.h, 'Color', 'white');
            orient landscape
        end
        
        function panel(o, index, varargin)
            p = inputParser;
            p.addRequired('index', @isscalar);
            p.addParameter('letter', char('a' + index - 1), @ischar);
            p.addParameter('xlab', '', @ischar);
            p.addParameter('ylab', '', @ischar);
            p.addParameter('title', '', @ischar);
            p.addParameter('x_shift', -0.08, @isscalar);
            p.addParameter('y_shift', 0.11, @isscalar);
            p.parse(index, varargin{:});
            r = p.Results;
            
            figure(o.h);
            ax = subplotp(o.Q, index);
            o.subplots{index} = ax;
            
            
            lettering(r.letter, r.x_shift, r.y_shift);
            
            o.xlabs{index} = r.xlab;
            xlabel(ax, r.xlab);
            o.ylabs{index} = r.ylab;
            ylabel(ax, r.ylab);
            o.titles{index} = r.title;
            title(ax, r.title);
        end
        
        function format(o)
            for index = 1:numel(o.subplots)
                ax = o.subplots{index};
                panel_format(ax);
                xlabel(ax, o.xlabs{index}, 'Interpreter', 'tex');
                ylabel(ax, o.ylabs{index}, 'Interpreter', 'tex');
                title(ax, o.titles{index}, 'Interpreter', 'tex');
            end
        end
        
        function print(o, dest, fname)
            figure(o.h);
            Utils.printto(dest, fname);
        end
    end
    
end