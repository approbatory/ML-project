function extractor_part(M, config, codename, idx_partition)
% Wrapper for EXTRACT for processing large movies

% Get start time
start_time = posixtime(datetime);
io_time = 0;

ABS_TOL = 1e-6;
SIGNAL_LOWER_THRESHOLD = 1e-6;
PARTITION_SIDE_LEN = 250;

if ~exist('config', 'var') || ~isfield(config, 'avg_cell_radius')
    error('"config.avg_cell_radius" must be specified.');
end

% Update config with defaults
config = get_defaults(config);

do_auto_partition=1;
if isfield(config, 'num_partitions_x') && ...
        isfield(config, 'num_partitions_y')
    do_auto_partition=0;
end
partition_overlap = ceil(config.avg_cell_radius * 2);

num_workers = 0;

% Delete existing parpool if using multi-gpu
if config.multi_gpu
    p = gcp('nocreate');
    delete(p);
end

% Override the gpu flag if necessary + handle multi-gpu case
if config.use_gpu && ~config.use_default_gpu
    dispfun(sprintf('%s: Getting GPU information... \n', datestr(now)),...
        config.verbose ~= 0);
    max_mem = 0;
    min_mem = inf;
    idx_max_mem = 0;
    c = gpuDeviceCount;
    gpuDevice([]);
    for idx_gpu = 1:c
        d = gpuDevice(idx_gpu);
        mem = d.AvailableMemory;
        dispfun(sprintf(...
            '\t \t \t GPU Device %d - %s: Available Memory: %.1f Gb\n', ...
            idx_gpu, d.Name, mem / 2^30), config.verbose ~= 0);
        min_mem = min(mem, min_mem);
        if mem > max_mem
            max_mem = mem;
            idx_max_mem = idx_gpu;
        end
    end
    if config.multi_gpu && c > 1
        avail_mem = min_mem;
        num_workers = c;
        % De-select last selected GPU
        gpuDevice([]);
        parpool('local', num_workers);
    else
        avail_mem = max_mem;
        gpuDevice(idx_max_mem);
        dispfun(sprintf('\t \t \t - Selecting GPU device %d \n', ...
            idx_max_mem), config.verbose ~= 0);
    end
    if c == 0
        warning('No GPU device was detected -- Setting use_gpu = 0 ');
        config.use_gpu = 0;
    end
end

% Set num_workers for parallel computation on CPUs
if ~config.use_gpu && config.parallel_cpu == 1
    % Default # of parallel workers is # cores -1
    num_workers = feature('numCores') - 1;
    if isfield(config, 'num_parallel_cpu_workers')
        if config.num_parallel_cpu_workers > num_workers + 1
            warning(['More parallel CPU workers than # of available ', ...
                'cores requested -- Using max available = %d'], ...
                num_workers + 1);
        else
            num_workers = config.num_parallel_cpu_workers;
            dispfun(sprintf(...
                '\t \t \t Setting up a pool with %d CPU workers \n', ...
                num_workers), config.verbose ~= 0);
        end
    end
    % Delete existing parpool if specs are different
    p = gcp('nocreate');
    if ~isempty(p)
        if p.NumWorkers ~=num_workers
            delete(p);
            parpool('local', num_workers);
        end
    else
        parpool('local', num_workers);    
    end
end

% Prevent plots when in parallel mode
if num_workers > 1
    config.plot_loss = 0;
end

[h, w, ~] = get_movie_size(M);

npt = config.num_frames;
% Determine the movie partitions
if ~do_auto_partition
    npx = config.num_partitions_x;
    npy = config.num_partitions_y;
else
    % Account for downsampling
    dss = config.downsample_space_by;
    if strcmp(dss, 'auto')
        dss = max(round(config.avg_cell_radius / ...
            config.min_radius_after_downsampling), 1);
    end
    % If estimation of full S is desired after downsampling, set dss = 1
    if config.reestimate_S_if_downsampled
        dss = 1;
    end
    config.downsample_space_by = dss;
    h_adjusted = h / dss;
    w_adjusted = w / dss;
    npx = max(round(w_adjusted / PARTITION_SIDE_LEN), 1);
    npy = max(round(h_adjusted / PARTITION_SIDE_LEN), 1);
    dispfun(sprintf(...
        '%s: Signal extraction will run on %d partitions (%dx%d) \n', ...
        datestr(now), npx * npy, npy, npx), config.verbose ~= 0);
end

% Get a circular mask (for movies with GRIN)
if config.crop_circular
    if ischar(M)
        error('To use the circular cropping feature, load the movie onto memory before calling EXTRACT.');
    else
        circular_mask = get_circular_mask(M);
        if isempty(config.movie_mask)
        	config.movie_mask = circular_mask;
        else
            % Apply user mask AND circular mask
            config.movie_mask = circular_mask & config.movie_mask;
        end
    end
end


num_partitions = npx * npy;
fov_occupation_total = zeros(h, w);
summary_image = zeros(h,w);
max_image = zeros(h,w);

summary = {};
S = {};
T = {};


% Divide movie into blocks and run EXTRACT
% parfor (idx_partition = 1:num_partitions, num_workers)
%%for idx_partition = num_partitions:-1:1
    dispfun(sprintf('%s: Signal extraction on partition %d (of %d):\n', ...
        datestr(now), idx_partition, num_partitions), config.verbose ~= 0);
    
    tic;
    % Get current movie partition from full movie
    [M_small, fov_occupation] = get_current_partition(...
        M, npx, npy, npt, partition_overlap, idx_partition);
    io_time = io_time + toc;

    % Sometimes partitions contain no signal. Terminate in that case
    std_M = nanstd(M_small(:));
    if (std_M < SIGNAL_LOWER_THRESHOLD) || (numel(M_small) == 0)
        dispfun('\t \t \t No signal detected, terminating...\n', ...
            config.verbose ==2);
	clear M_small;
	save(sprintf('extract_part_for_%s_index_%d.mat', codename, idx_partition));
	return;
    end
    config_this = config;
    % If S_init is given, feed only part of it consistent with partition
    if ~isempty(config_this.S_init)
        S_init = config.S_init(fov_occupation(:), :);
        S_init(:, sum(S_init, 1)<=ABS_TOL) = [];
        config_this.S_init = S_init;
    end
    % Distribute mask to partitions
    if ~isempty(config_this.movie_mask)
        [h_this, w_this, ~] = size(M_small);
        config_this.movie_mask = config_this.movie_mask(fov_occupation(:));
        config_this.movie_mask = reshape(config_this.movie_mask, h_this, w_this);
    end
    % Run EXTRACT for current partition
    [S_this, T_this, summary_this] = run_extract(M_small, config_this);
    dispfun(sprintf('\t \t \t Count: %d cells.\n', ...
        size(S_this, 2)), config.verbose ~= 0);
    S_this = sparse(double(S_this));
    % Un-trim the pixels
    S_temp = spalloc(h * w, size(S_this, 2), nnz(S_this));
    S_temp(fov_occupation(:), :) = S_this;
    S_this = S_temp;

    % Update summary image & max_image
    if isfield(summary_this, 'summary_image')
        summary_image(fov_occupation(:)) = summary_this.summary_image;
        max_image(fov_occupation(:)) = summary_this.max_image;
    else
        summary_image(fov_occupation(:)) = max(M_small, [], 3);
        max_image(fov_occupation(:)) = max(M_small, [], 3);
    end
    summary_this.fov_occupation = fov_occupation;
    summary{idx_partition} = summary_this;
    if ~isempty(S_this)
        S{idx_partition} = S_this;
        T{idx_partition} = T_this';
    end
    fov_occupation_total = fov_occupation_total + fov_occupation;
%%end

clear M_small S T S_temp
%S_this = sparse(double(S_this));
save(sprintf('extract_part_for_%s_index_%d.mat', codename, idx_partition));
return;

% Concatenate components from different partitions
S = cell2mat(S(~cellfun(@isempty, S)));
T = cell2mat(T(~cellfun(@isempty, T)));
summary = [summary{~cellfun(@isempty, summary)}];

dispfun(sprintf('%s: Total of %d cells are found.\n', ...
    datestr(now),size(S,2)),config.verbose~=0);

if config.remove_duplicate_cells
    dispfun(sprintf('%s: Removing duplicate cells...\n', ...
    datestr(now), size(S, 2)), config.verbose ~= 0);

    overlap_idx = find(fov_occupation_total - 1);
    if ~isempty(S)
        idx_trash = find_duplicate_cells(S, T, overlap_idx);
        S(:, idx_trash) = [];
        T(:, idx_trash) = [];
    end

    dispfun(sprintf(...
        '%s: %d cells were retained after removing duplicates.\n', ...
        datestr(now), size(S, 2)), config.verbose ~=0);
end

% Get total runtime (minus time spent reading from disk)
end_time = posixtime(datetime);
total_runtime = end_time - start_time - io_time;

info.version = '0.6.1';
info.summary = summary;
info.runtime = total_runtime;
info.summary_image = summary_image;
info.max_image = max_image;

output.spatial_weights = reshape(S, h, w, size(S, 2));
output.temporal_weights = T;
output.info = info;
output.config = config;

fprintf('%s: All done! \n', datestr(now));

end
