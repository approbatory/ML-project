function output = extractor_aggregator(varargin)
[S, T, summary] = deal(cell(1,nargin));

first_one = true;
for i = 1:nargin
    fname = varargin{i};
    L = load(fname);
    %assert(isfield(L, 'S_this'), '%s does not contain S_this', fname);
    %assert(isfield(L, 'T_this'), '%s does not contain T_this', fname);
    %assert(isfield(L, 'summary_this'), '%s does not contain summary_this', fname);
    if ~isfield(L, 'S_this') || ~isfield(L, 'T_this') || ~isfield(L, 'summary_this')
        fprintf('Skipping file %s\n', fname);
        continue;
    end
    S{i} = L.S_this;
    T{i} = L.T_this';
    summary{i} = L.summary_this;
    
    if first_one
        fov_occupation_total = L.fov_occupation;
        summary_image = L.summary_image;
        max_image = L.max_image;
        config = L.config;
        w = L.w;
        h = L.h;
        first_one = false;
    else
        fov_occupation_total = fov_occupation_total + L.fov_occupation;
        summary_image = summary_image + L.summary_image;
        max_image = max(max_image, L.max_image);
    end     
end

% Concatenate components from different partitions
S = cell2mat(S(~cellfun(@isempty, S)));
T = cell2mat(T(~cellfun(@isempty, T)));
summary = [summary{~cellfun(@isempty, summary)}];


if config.remove_duplicate_cells
    dispfun(sprintf('%s: Removing duplicate cells...\n', ...
    datestr(now), size(S, 2)), config.verbose ~= 0);

    overlap_idx = find(fov_occupation_total == 2);
    if ~isempty(S)
        idx_trash = find_duplicate_cells_sparse(S, T, overlap_idx);
        S(:, idx_trash) = [];
        T(:, idx_trash) = [];
    end

    dispfun(sprintf(...
        '%s: %d cells were retained after removing duplicates.\n', ...
        datestr(now), size(S, 2)), config.verbose ~=0);
end


info.version = '0.6.1';
info.summary = summary;
info.summary_image = summary_image;
info.max_image = max_image;

output.spatial_weights = reshape(ndSparse(S), h, w, size(S,2));%ndSparse(S, h, w, size(S,2));
output.temporal_weights = T;
output.info = info;
output.config = config;

fprintf('%s: All done! \n', datestr(now));
