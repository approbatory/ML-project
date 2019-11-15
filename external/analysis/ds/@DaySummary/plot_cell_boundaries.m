function plot_cell_boundaries(obj, varargin)
    % Default display options
    show_cell_map_background = 1;
    cells_only = 0;
    linespec = [];
    linewidth = 1;
    for k = 1:length(varargin)
        if ischar(varargin{k})
            switch lower(varargin{k})
                case 'nobackground'
                    show_cell_map_background = 0;
                case {'cells', 'cellsonly'}
                    cells_only = 1;
                case 'linespec'
                    % Overrides default coloring by cell label
                    linespec = varargin{k+1};
                case 'linewidth'
                    linewidth = varargin{k+1};
            end
        end
    end

    % Display the cell map
    if show_cell_map_background
        imagesc(obj.cell_map_ref_img);
        colormap gray;
    else
        [h, w] = size(obj.cell_map_ref_img);
        plot([1 w w 1], [h h 1 1], 'k');
        set(gca, 'YDir', 'Reverse');
    end
    hold on;
    axis equal tight;
            
    for k = 1:obj.num_cells
        boundary = obj.cells(k).boundary;
        
        % By default, boundary colors are defined by cell label, but can be
        % overriden by 'linespec' varargin
        if isempty(linespec)
            if obj.is_cell(k)
                ls = 'g';
            else
                ls = 'r';
            end
        else
            ls = linespec;
        end
        
        if obj.is_cell(k) || ~cells_only
            plot(boundary(:,1), boundary(:,2), ls, 'linewidth', linewidth);
        end
        
    end
    hold off;
end