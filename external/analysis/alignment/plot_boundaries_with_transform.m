function plot_boundaries_with_transform(ds, linespec, linewidth, filled_cells, tform)
    % Plot boundaries as a single color, with an optional transform. Can
    % subselect cells to be filled in
    if (nargin < 3)
        linewidth = 1;
        filled_cells = [];
        tform = [];
    elseif (nargin < 4)
        filled_cells = [];
        tform = [];
    elseif (nargin < 5)
        tform = [];
    end
    
    if isempty(linewidth)
        linewidth = 1;
    end
    
    for k = 1:ds.num_cells
        boundary = ds.cells(k).boundary;
        if ~isempty(tform) % Optional spatial transform
            boundary = transformPointsForward(tform, boundary);
        end
        if ismember(k, filled_cells)
            fill(boundary(:,1), boundary(:,2), linespec,...
                 'LineWidth', linewidth,...
                 'FaceAlpha', 0.4);
        elseif ds.is_cell(k)
            plot(boundary(:,1), boundary(:,2), linespec, 'LineWidth', linewidth, 'HitTest', 'off');
        end
        hold on;
    end
    axis equal tight;
    set(gca, 'YDir', 'Reverse');
end