function browse_rasters_touch(ds, varargin)
% Tool for browsing single cell rasters of a single day (i.e. DaySummary),
% but without keyboard interaction!

h_fig = [];
page_idx = [];
cell_idx = [];
for i = 1:length(varargin)
    vararg = varargin{i};
    if ischar(vararg)
        switch lower(vararg)
            case 'fig'
                h_fig = varargin{i+1};
            case {'page', 'page_idx'}
                page_idx = varargin{i+1};
            case {'cell', 'cell_idx'}
                cell_idx = varargin{i+1};
        end
    end
end

if isempty(h_fig)
    h_fig = figure;
else
    figure(h_fig);
    clf;
end

% Main drawing routine
%------------------------------------------------------------
cells_per_page = [2 4];
num_cells_per_page = prod(cells_per_page);

cell_indices = find(ds.is_cell);
cell2page = ceil((1:length(cell_indices))/num_cells_per_page);
max_page = max(cell2page);

if isempty(page_idx)
    if isempty(cell_idx)
        page_idx = 1;
    else
        page_idx = cell2page(cell_indices==cell_idx);
    end
end

page_idx = max(1, page_idx);
page_idx = min(max_page, page_idx);

cells_on_page = cell_indices(cell2page==page_idx);

for i = 1:length(cells_on_page)
    cell_idx = cells_on_page(i);
    
    h_sp = subplot(cells_per_page(1), cells_per_page(2), i);
    ds.plot_cell_raster(cell_idx);
    title(sprintf('Cell %d', cell_idx));
    
    set(h_sp, 'ButtonDownFcn', {@raster_clicked, cell_idx});
end

% Navigation controls
%------------------------------------------------------------
prev_btn = uicontrol('Style', 'pushbutton',...
    'String', '<<',...
    'Units', 'normalized',...
    'Position', [0 0.1 0.05 0.9],...
    'Callback', {@page_button_clicked, page_idx-1});
if page_idx == 1
    prev_btn.Enable = 'off';
end

next_btn = uicontrol('Style', 'pushbutton',...
    'String', '>>',...
    'Units', 'normalized',...
    'Position', [0.95 0 0.05 1],...
    'Callback', {@page_button_clicked, page_idx+1});
if page_idx == max_page
    next_btn.Enable = 'off';
end

jump_btn = uicontrol('Style', 'pushbutton',...
    'String', 'C#',...
    'Units', 'normalized',...
    'Position', [0 0 0.05 0.1],...
    'Callback', @jump_to_cell);

    function jump_to_cell(~, ~)
        prompt = sprintf('Enter cell index [1-%d]', ds.num_cells);
        jump_idx = inputdlg(prompt, 'Jump to cell');
        jump_idx = str2double(jump_idx{1});
        if ~isempty(jump_idx)
            if (1 <= jump_idx) && (jump_idx <= ds.num_cells)
%                 browse_rasters_touch(ds, 'fig', h_fig, 'cell_idx', jump_idx);
                view_raster_touch(ds, jump_idx, 'fig', h_fig);
            end
        end
    end

    function page_button_clicked(~, ~, page_idx)
        browse_rasters_touch(ds, 'fig', h_fig, 'page_idx', page_idx);
    end

    function raster_clicked(~, ~, cell_idx)
        view_raster_touch(ds, cell_idx, 'fig', h_fig);
    end

end % draw_raster_page
