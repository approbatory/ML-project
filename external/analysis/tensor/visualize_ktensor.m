function [Ax, G, BigAx] = visualize_ktensor(X, varargin)
% PLOT Plot the factor matrices of a ktensor
%
%     [Ax, G, BigAx] = PLOT(X)
%
%MATLAB Tensor Toolbox.
%Copyright 2015, Sandia Corporation.

% This is the MATLAB Tensor Toolbox by T. Kolda, B. Bader, and others.
% http://www.sandia.gov/~tgkolda/TensorToolbox.
% Copyright (2015) Sandia Corporation. Under the terms of Contract
% DE-AC04-94AL85000, there is a non-exclusive license for use of this
% work by or on behalf of the U.S. Government. Export of this data may
% require a license from the United States Government.
% The full license terms can be found in the file LICENSE.txt


% check inputs
if ~isa(X,'ktensor')
    error('input must be a ktensor')
end

% tensor order, shape, and rank
nf = ndims(X);
sz = size(X);
nr = length(X.lambda);

% parse optional inputs
params = inputParser;

params.addParameter('nfactors', nr);
params.addParameter('align', true);
params.addParameter('space', 0.1);
params.addParameter('title', cell(1,nf));
params.addParameter('plots', repmat({'line'}, [1 nf]));
params.addParameter('a', repmat({10}, [1 nf]));
params.addParameter('c', cell(1,nf));
params.addParameter('markertype', repmat({'o'}, [1 nf]));
params.addParameter('filled', true(1,nf));
params.addParameter('linespec', repmat({'-'}, [1 nf]));
params.addParameter('link_yax', false(1,nf));
params.addParameter('ylims', cell(1, nf));
params.addParameter('greedy', nr>7);
params.addParameter('linewidth', ones(1,3));
params.addParameter('xlabel', repmat({''}, [1 nf]));
params.addParameter('XTick', cell(1,nf));
params.addParameter('XTickLabel', cell(1,nf));
params.addParameter('axes', []);

params.parse(varargin{:});
res = params.Results;

% set up the axes
if isempty(res.axes)
    [Ax,BigAx] = setup_axes(nr, nf, res.space);
else
    Ax = res.axes;
    if size(Ax,1) ~= nr
        error('User-provided Axes do not match the rank of the provided tensor.')
    elseif size(Ax,2) ~= nf
        error('User-provided Axes do not match the order of the provided tensor.')
    end
    BigAx = []; % don't return if user provides axes
end
format_axes(Ax, res.title, res.xlabel, res.XTick, res.XTickLabel);

% main loop %
% iterate over factors (columns of plot)
G = gobjects(nr, nf);
for f = 1:nf

    % iterate over model rank (rows of plot)
    for r = 1:nr
        
        % fetch the axes to plot on
        axes(Ax(r,f));
        hold on

        % plot x vs y
        x = 1:sz(f);
        y = X.u{f}(:, r);

        % make the plot
        switch res.plots{f}
            case 'line'
                G(r, f) = plot(x, y, res.linespec{f}, 'linewidth', res.linewidth(f));

            case 'scatter'
                if isempty(res.c{f})
                    if res.filled
                        G(r, f) = scatter(x, y, res.a{f}, res.markertype{f}, 'filled');
                    else
                        G(r, f) = scatter(x, y, res.a{f}, res.markertype{f});
                    end
                else
                    if res.filled
                        G(r, f) = scatter(x, y, res.a{f}, res.c{f}, res.markertype{f}, 'filled');
                    else
                        G(r, f) = scatter(x, y, res.a{f}, res.c{f}, res.markertype{f});
                    end
                end

            case 'bar'
                G(r, f) = bar(x, y);

            otherwise
                error('Did not understand plot type.')

        end

        % make axes tight
        axis tight
    end
end


% make the ylims look nice before returning
pretty_ylims(Ax, res.ylims, res.link_yax);        

%%%%%%%%%%%%%%%%%%%
% LOCAL FUNCTIONS %
%%%%%%%%%%%%%%%%%%%

function [Ax,BigAx] = setup_axes(nr, nf, space)

    % allocate storage
    Ax = gobjects(nr,nf);
    BigAx = gobjects(1,nf);
    
    % setup axes
    for f = 1:nf
        % invisible subplot bounding box
        BigAx(f) = subplot(1, nf, f);
        set(BigAx(f),'Visible','off')
        pos = get(BigAx(f),'Position');
        w = pos(3);
        h = pos(4)/nr;
        pos(1:2) = pos(1:2) + space*[w h];

        % subaxes
        for r = 1:nr
            axPos = [pos(1) pos(2)+(nr-r)*h w*(1-space) h*(1-space)];
            Ax(r,f) = axes('Position',axPos);
        end
    end

function format_axes(Ax, ttl, xlab, xt, xtl)
    [nr,nf] = size(Ax);
    for f = 1:nf
        for r = 1:nr
            if ~isempty(ttl{f}) && r == 1
                set(Ax(1,f).Title,'String',ttl{f})
            end
            if r ~= nr
                set(Ax(r,f),'XTick',[]);
                set(Ax(r,f),'XColor','w');
            else
                if ~isempty(xt{f})
                    Ax(r,f).XTick = xt{f};
                end
                if ~isempty(xtl{f})
                    Ax(r,f).XTickLabel = xtl{f};
                end
                xlabel(Ax(r,f),xlab{f})
            end
            if mod(r,2) == 0
                set(Ax(r,f),'YAxisLocation','right')
            end
        end
    end

function pretty_ylims(Ax, ylimits, link)

    % dimensions
    [nr,nf] = size(Ax);

    % set yticks %
    for f = 1:nf
        % if ylim pre-specified for this column
        if ~isempty(ylimits{f})
            yl = ylimits{f};
            yt = pretty_axticks(yl);
            set(Ax(:,f), 'ylim', yl, 'YTick', yt);

        % if ylim not pre-specified but linked in this column
        elseif link(f)
            % find new ylims
            yl = [Inf -Inf];
            for r = 1:nr
                ylr = get(Ax(r,f), 'ylim');
                if yl(1) > ylr(1)
                    yl(1) = ylr(1);
                end
                if yl(2) < ylr(2)
                    yl(2) = ylr(2);
                end
            end

            % set ylim and yticks
            yt = pretty_axticks(yl);
            set(Ax(:,f), 'ylim', yl, 'YTick', yt);

        % if ylim not pre-specified, just make the ticks pretty
        else
            for r = 1:nr
                yt = pretty_axticks(get(Ax(r,f), 'ylim'));
                set(Ax(r,f), 'YTick', yt);
                set(Ax(r,f), 'TickDir', 'out')
            end
        end
    end

function ryl = pretty_axticks(yl)
    t0 = round(ceil(yl(1)*100)/100,2);
    t1 = round(floor(yl(2)*100)/100,2);
    ryl = [t0 t1];
