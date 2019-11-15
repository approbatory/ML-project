function browse_rasters(ds, varargin)
% Tool for browsing single cell rasters of a single day (i.e. DaySummary)
%
% Usage:
%     m3d15 = DaySummary(sources.maze, 'rec001', 'reconst', 'excludeprobe');
%       14-Aug-2015 13:58:23: Loaded data from rec001/rec_150813-115003.mat
%       14-Aug-2015 13:58:23: Loaded classification from rec001/class_150813-115320.txt
%     browse_rasters(m3d15);
%

cell_indices = [];

for k = 1:length(varargin)
    if ischar(varargin{k})
        switch lower(varargin{k})
            case 'cells' % Specify the IDs of cells to browse through
                cell_indices = varargin{k+1};
        end
    end
end

% By default, browse through rasters of all classified cells
if isempty(cell_indices)
    cell_indices = find(ds.is_cell);
end
num_cells = length(cell_indices);

% Display settings
cells_per_page = [2 4];
num_cells_per_page = prod(cells_per_page);
num_pages = ceil(num_cells / num_cells_per_page);

h = figure;

page_idx = 1;
while (1)
    draw_page(page_idx);

    % Ask user for command
    prompt = sprintf('Raster browser (page %d of %d) >> ', page_idx, num_pages);
    resp = strtrim(input(prompt, 's'));
    
    val = str2double(resp);
    if (~isnan(val)) % Is a number (cell index)
        if ismember(val, cell_indices)
            fprintf('  Cell %d selected!\n', val);
            view_detailed_raster(ds, val);
            
            % When done, return to the page that contains the cell
            page_idx = ceil(find(cell_indices==val,1) / num_cells_per_page);
        else
            fprintf('  Sorry, %d is not a valid cell index\n', val);
        end
    else
        resp = lower(resp);
        switch (resp)    
            case {'n', ''} % Next page
                if (page_idx < num_pages)
                    page_idx = page_idx + 1;
                else
                    fprintf('  Already at final page!\n');
                end

            case 'p' % Previous page
                if (page_idx > 1)
                    page_idx = page_idx - 1;
                else
                    fprintf('  Already at first page!\n');
                end

                
            case 'f' % First page
                page_idx = 1;
                
            case 'q' % Exit
                close(h);
                break;

            otherwise
                fprintf('  Could not parse "%s"\n', resp);

        end % switch
    end
end

    % Helper functions
    %------------------------------------------------------------
    function draw_page(page_idx)
        clf;
        
        cells_on_page = get_cells_on_page(page_idx);
        for i = 1:length(cells_on_page)
            subplot(cells_per_page(1), cells_per_page(2), i);
        
            cell_idx = cells_on_page(i);
            ds.plot_cell_raster(cell_idx);
            title(sprintf('Cell %d', cell_idx));
        end
        
        function cells_on_page = get_cells_on_page(page_idx)
            index_to_cells = [1+num_cells_per_page*(page_idx-1) num_cells_per_page*page_idx];
            index_to_cells(2) = min(index_to_cells(2), num_cells);
            cells_on_page = cell_indices(index_to_cells(1):index_to_cells(2));
        end % get_cells_on_page
        
    end % draw_page

end % browse_rasters