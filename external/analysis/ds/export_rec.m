function export_rec(ds, varargin)
% Export the contents of DaySummary into a new REC file, with specified
% manipulations to the underlying data. Only classified cells will be
% exported.
%
% NOTE: Do not export from a DaySummary that has been initialized with
% 'excludeprobe' as this can lead to synchronization problems with the
% original trial metadata!

% Cell selection options. By default, export only classified cells
cell_indices = find(ds.is_cell);

% Filter options
truncate_filter = false;

% Trace options
remove_baseline = false;
normalize_trace = false;

for k = 1:length(varargin)
    vararg = varargin{k};
    if ischar(vararg)
        switch lower(vararg)
            case {'fix', 'fix_baseline', 'remove_baseline'}
                fprintf('  Traces will be baseline corrected...\n');
                remove_baseline = true;
                
            case 'norm'
                fprintf('  Traces will be normalized...\n');
                normalize_trace = true;
                
            case 'truncate'
                fprintf('  Filters will be truncated...\n');
                truncate_filter = true;
                
            case {'keep_unlabeled', 'remove_noncells'}
                fprintf('  Removing cells classified to be not a cell (but keeping unlabeled ones)...\n');
                labels = {ds.cells.label};
                not_a_cell = cellfun(@strcmp,...
                    labels, repmat({'not a cell'}, size(labels)));
                cell_indices = find(~not_a_cell);
        end
    end
end

num_cells = length(cell_indices);
num_frames = ds.full_num_frames;

% Confirm that probe trials have not been excluded, which can lead to
% synchronization problems for the resulting REC file
trace_len = length(ds.get_trace(1));
assert(num_frames == trace_len,...
       'Error! Length of retrieved traces (%d) does not match full_num_frames (%d)!\n',...
       trace_len, num_frames);
traces = zeros(num_frames, num_cells);
   
[height, width] = size(ds.cells(1).im);
filters = zeros(height, width, num_cells);

for k = 1:num_cells
    cell_idx = cell_indices(k);
    filters(:,:,k) = ds.cells(cell_idx).im;
    if truncate_filter
        filters(:,:,k) = ds.cells(cell_idx).mask .* filters(:,:,k);
    end
    
    trace = ds.get_trace(cell_idx);
    if remove_baseline
        trace = fix_baseline(trace);
    end
    if normalize_trace
        tscale = [min(trace) max(trace)];
        trace = (trace - tscale(1))/diff(tscale);
    end
    traces(:,k) = trace;
end

% Save as Rec file
%------------------------------------------------------------
info.type = 'export_rec';
info.num_pairs = num_cells;

% Note the parameters used in export
info.options.truncate_filter = truncate_filter;
info.options.remove_baseline = remove_baseline;
info.options.normalize_trace = normalize_trace;

timestamp = datestr(now, 'yymmdd-HHMMSS');
rec_savename = sprintf('rec_%s.mat', timestamp);
save(rec_savename, 'info', 'filters', 'traces', '-v7.3');