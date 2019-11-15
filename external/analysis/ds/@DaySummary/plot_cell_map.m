function plot_cell_map(obj, color_grouping, varargin)
    % Optional argument allows for specification of color used for
    % the cell in the cell map. The color specification is defined
    % as follows:
    %   color_grouping = {[1, 2, 3, 4], 'w';
    %                     [5, 6], 'r';
    %                     [10], 'g'}
    % means that cells [1, 2, 3, 4] will be displayed in white,
    % cells [5, 6] in red, and [10] in green.

    enable_class_colors = 0;
    for k = 1:length(varargin)
        vararg = varargin{k};
        if ischar(vararg)
            switch lower(vararg)
                % By default, if 'color_grouping' is provided, the
                % default boundary coloring by classification is
                % disabled. The following toggle enables the
                % class colors, which is then overridden by the
                % specified 'color_grouping'
                case 'enable_class_colors'
                    enable_class_colors = 1;
            end
        end
    end

    cell_colors = cell(obj.num_cells, 1);

    % By default, color the cells based on classification
    if ~exist('color_grouping', 'var') || enable_class_colors
        for k = 1:obj.num_cells
            % Note: Unlabeled cells remain white!
            if isempty(obj.cells(k).label)
                cell_colors{k} = 'w';
            else
                if obj.is_cell(k)
                    cell_colors{k} = 'g';
                else
                    cell_colors{k} = 'r';
                end
            end
        end
    end

    % Color grouping provided, unpack into cell_linespec
    if exist('color_grouping', 'var')
        num_groups = size(color_grouping, 1);
        for group = 1:num_groups
            cells_in_group = color_grouping{group,1};
            if iscolumn(cells_in_group)
                cells_in_group = cells_in_group'; % Want row
            end
            for cell_idx = cells_in_group
                cell_colors{cell_idx} = color_grouping{group,2};
            end
        end
    end

    % Display the cell map
    imagesc(obj.cell_map_ref_img, 'HitTest', 'off');
    colormap gray;
    axis equal tight;

    hold on;
    for k = 1:obj.num_cells
        if ~isempty(cell_colors{k})
            boundary = obj.cells(k).boundary;

            if strcmp(cell_colors{k}, 'w')
                linewidth = 0.25;
                alpha = 0.1; % Lower is more transparent
            else
                linewidth = 1.0;
                alpha = 0.3;
            end

            % Note: We are embedding the cell index in the
            %   z-dimension of the graphics object (hack!)
            fill3(boundary(:,1), boundary(:,2),...
                 k*ones(size(boundary,1),1),...
                 cell_colors{k},...
                 'EdgeColor', cell_colors{k},...
                 'LineWidth', linewidth,...
                 'FaceAlpha', alpha);
        end
    end
    hold off;

    dcm = datacursormode(gcf);
    set(dcm, 'DisplayStyle', 'datatip');
    set(dcm, 'UpdateFcn', @cellmap_tooltip);
    datacursormode on;
end