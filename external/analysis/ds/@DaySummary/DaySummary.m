% Summary of PlusMaze data for a single day.
%
% Inputs:
%   ds_source: Three possibilities:
%       1) Plus maze text file
%       2) Struct 's' containing the following fields:
%           - s.maze: Path to plus maze text file (required)
%           - s.behavior: Path to behavioral video (optional; e.g. mp4)
%           - s.tracking: Path to tracking text file (optional; *.xy)
%       3) Empty string. Fake trial metadata will be generated that matches
%           the number of frames in the provided rec file. A hack to allow
%           use of the 'analysis' codebase with non-PlusMaze datasets.
%
%   rec_dir: Directory containing 
%       - Filters and traces in a "rec_*.mat" file (required)
%       - Classification results in a "class_*.txt" file (optional)
%
% Output:
%   DaySummary object
%
% Example usage:
%   ds = DaySummary('c11m1d12_ti2.txt', 'rec001');
%
classdef DaySummary < handle
    properties
        cells
        trials
        switchdata
        
        num_cells
        num_trials
        
        trial_indices
        full_num_frames
    end
    
    properties (SetAccess = private, Hidden=true)
        trace_corrs
        trace_range
        cell_distances
        cell_map_ref_img
        behavior_vid
        behavior_ref_img
        is_tracking_loaded
        is_eventdata_loaded
        orig_trial_indices
    end
        
    methods
        function obj = DaySummary(ds_source, rec_dir, varargin)
            % Handle optional input
            exclude_probe_trials = 0;
            for k = 1:length(varargin)
                if ischar(varargin{k})
                    switch lower(varargin{k})
                        case {'excludeprobe', 'noprobe'}
                            exclude_probe_trials = 1;
                    end
                end
            end
            
            % Load extraction data (i.e. filters & traces)
            %------------------------------------------------------------
            data_source = get_most_recent_file(rec_dir, 'rec_*.mat');
            data = load(data_source);
            obj.num_cells = data.info.num_pairs;
            fprintf('%s: Loaded %d filters and traces (%s) from %s\n',...
                    datestr(now), obj.num_cells, data.info.type, data_source);
                
            trace_num_frames = size(data.traces, 1);
            
            % Check if DaySummary session data is provided as a struct
            %------------------------------------------------------------
            if isstruct(ds_source)
                plusmaze_txt = ds_source.maze;
            else
                plusmaze_txt = ds_source;
            end
                       
            if ~isempty(plusmaze_txt)
                % PlusMaze metadata provided. Read trial metadata
                [trial_indices, loc_info, trial_durations] = parse_plusmaze(plusmaze_txt); %#ok<*PROP>
                fprintf('%s: Loaded trial metadata from %s\n', datestr(now), plusmaze_txt);
            else
                % PlusMaze metadata NOT provided. Make up fake info so that
                % DaySummary object can be instantiated anyway.
                trial_indices = [1 2 3 trace_num_frames];
                loc_info = {'east', 'north', 'north'};
                trial_durations = 10.0;
                fprintf('%s: Generated fake trial metadata to match number of frames (%d) in rec file (%s)!\n',...
                    datestr(now), trace_num_frames, data_source);
            end
            
            % Check that the length of traces is consistent with the table
            % of trial indices.
            obj.full_num_frames = trial_indices(end,end);
            assert(trace_num_frames == obj.full_num_frames,...
                'Error: Length of traces does not match trial index table!');
            
            % Optional exclusion of probe trials. Effectively, we are
            % "deleting" the lines of the plus maze text file that
            % correspond to probe trials.
            if (exclude_probe_trials)
                is_probe = strcmp(loc_info(:,1), 'north') | ...
                           strcmp(loc_info(:,1), 'south');
                       
                trial_indices = trial_indices(~is_probe,:);
                loc_info = loc_info(~is_probe,:);
                trial_durations = trial_durations(~is_probe);
            end
            
            % Parse by TRIAL
            %------------------------------------------------------------
            num_trials = size(trial_indices, 1);
            turns = cell(num_trials, 1);
            traces = cell(num_trials, 1);
            centroids = cell(num_trials, 1);

            for k = 1:num_trials
                trial_frames = trial_indices(k,1):...
                               trial_indices(k,end);

                num_frames_in_trial = length(trial_frames);
                traces{k} = data.traces(trial_frames, :)';
                turns{k} = obj.compute_turn(loc_info{k,1}, loc_info{k,3});
                centroids{k} = zeros(num_frames_in_trial, 2);
            end
            
            obj.num_trials = num_trials;
            
            % NOTE: We have to be quite careful about how we handle trial
            % indices if we apply probe trial elimination.
            %
            % Basically, when we need to refer to the raw data, we need the
            % _original_ frame indices, since the raw data sources do not
            % omit probe trials.
            %
            % On the other hand, when we later access content from the
            % probe-trial removed DaySummary (e.g. get_trace), the
            % resulting content does not know about the omitted trials.
            % Hence, for self-consistency of the DaySummary instance, we
            % need to "compress" the trial indices to remove gaps due to 
            % the omitted probe trials.
            %
            % Note that we _still_ need 'full_num_frames' in addition to
            % 'orig_trial_indices'. Namely, we can't always assume that
            % orig_trial_indices(end,end) == full_num_frames; this equality
            % is false when the session ends on a probe trial, and probe
            % trials are eliminated by DaySummary instantiation.
            %
            obj.trial_indices = compress_frame_indices(trial_indices, [0 0]);
            obj.orig_trial_indices = double(trial_indices);
            obj.trials = struct(...
                'start', loc_info(:,1),...
                'goal',  loc_info(:,2),...
                'end',   loc_info(:,3),...
                'correct', cellfun(@strcmp, loc_info(:,2), loc_info(:,3), 'UniformOutput', false),...
                'turn',  turns,...
                'time',  num2cell(trial_durations),...
                'traces', traces,...
                'events', [],...
                'centroids', centroids);
            
            full_traces = cell2mat(traces'); % [num_cells x num_frames]
            fprintf('  Computing auxiliary trace parameters...');
            tic;
            obj.trace_corrs = corr(full_traces');
            obj.trace_range = [min(full_traces,[],2) max(full_traces,[],2)];
            t = toc;
            fprintf(' Done (%.1f sec)\n', t);
            
            % Parse by CELL
            %------------------------------------------------------------
            class = cell(obj.num_cells,1);
            
            images = squeeze(num2cell(data.filters, [1 2])); % images{k} is the 2D image of cell k
            boundaries = cell(size(images));
            masks = cell(size(images));
            coms = cell(size(images)); % Center of mass

            fprintf('  Computing auxiliary spatial parameters...');
            tic;
            [height, width] = size(images{1});
            for k = 1:obj.num_cells
                boundary = compute_ic_boundary(images{k}, 0.3);
                boundaries{k} = boundary{1}; % Keep only the longest boundary!
                masks{k} = poly2mask(boundaries{k}(:,1), boundaries{k}(:,2), height, width);
                
                % Compute the center of mass
%                 masked_filter = masks{k}.*images{k};
                masked_filter = images{k}; % Masking sometimes yields numerical instabilities
                com = [(1:width)*sum(masked_filter,1)';
                       (1:height)*sum(masked_filter,2)];
                coms{k} = com / sum(masked_filter(:));
            end
            t = toc;
            fprintf(' Done (%.1f sec)\n', t);
            
            obj.cells = struct(...
                'im', images,...
                'boundary', boundaries,...
                'mask', masks,...
                'com', coms,...
                'label', class);
            
            % Compute distances among all sources
            fprintf('  Computing distances between all sources...');
            tic;
            D = Inf*ones(obj.num_cells);
            for i = 1:(obj.num_cells-1)
                for j = (i+1):obj.num_cells
                    delta = obj.cells(i).com - obj.cells(j).com;
                    D(i,j) = norm(delta);
                end
            end
            obj.cell_distances = min(D, D'); % Make symmetric. Note Infs.
            t = toc;
            fprintf(' Done (%.1f sec)\n', t);
            
            % Precompute cell map image, to avoid doing it each time
            [height, width] = size(obj.cells(1).im);
            ref_image = zeros(height, width);
            for k = 1:obj.num_cells
                ref_image = ref_image + obj.cells(k).im;
            end
            obj.cell_map_ref_img = ref_image;
            
            % Load classification, if available
            %------------------------------------------------------------
            class_source = get_most_recent_file(rec_dir, 'class_*.txt');
            if ~isempty(class_source)
                obj.load_class(class_source);
                fprintf('%s: Loaded classification from %s\n', datestr(now), class_source);
            end
                       
            % Other initialization
            %------------------------------------------------------------
            obj.behavior_vid = [];
            obj.is_tracking_loaded = false;
            
            % Load behavior video and tracking data if available
            if isstruct(ds_source)
                if isfield(ds_source, 'behavior')
                    obj.load_behavior_movie(ds_source.behavior);
                end
                if isfield(ds_source, 'tracking')
                    obj.load_tracking(ds_source.tracking);
                end
            end
            
            % Event data
            obj.is_eventdata_loaded = false;
            event_source = get_most_recent_file(rec_dir, 'events_*.mat');
            if ~isempty(event_source)
                obj.load_events(event_source);
            end
            
            % Fill in switch data
            obj.switchdata = struct(...
                'pre_switch_trials', [],...
                'post_switch_trials', [],...
                'changing_path_start', '',...
                'constant_path_start', '',...
                'changing_path_cutoff', []);
        end
        
        % Helper functions
        %------------------------------------------------------------        
        function turn = compute_turn(~, start, final)
            % TODO: Turn into Static
            path = {start, final};
            if (all(strcmp(path, {'east', 'south'})) || ...
                all(strcmp(path, {'south', 'west'})) || ...
                all(strcmp(path, {'west', 'north'})) || ...
                all(strcmp(path, {'north', 'east'})))
                turn = 'left';
            else
                turn = 'right';
            end
        end % compute_turn

        function filtered_trials = filter_trials(obj, varargin)
            filtered_trials = true(1, obj.num_trials);
            for k = 1:length(varargin)
                vararg = varargin{k};
                if ischar(vararg)
                    switch lower(vararg)
                        case {'range', 'inds', 'trial_inds'}
                            % Convert list of trial indices to logical vector
                            lv = ismember(1:obj.num_trials, varargin{k+1});
                            filtered_trials = filtered_trials & lv;
                        case 'incorrect'
                            filtered_trials = filtered_trials &...
                                (~strcmp({obj.trials.goal}, {obj.trials.end}));
                        case 'correct'
                            filtered_trials = filtered_trials &...
                                strcmp({obj.trials.goal}, {obj.trials.end});
                        case 'start'
                            filtered_trials = filtered_trials &...
                                trial_filter(obj, {obj.trials.start}, varargin{k+1});
                        case 'end'
                            filtered_trials = filtered_trials &...
                                trial_filter(obj, {obj.trials.end}, varargin{k+1});
                        case 'turn'
                            filtered_trials = filtered_trials &...
                                trial_filter(obj, {obj.trials.turn}, varargin{k+1});
                    end
                end
            end
            filtered_trials = filtered_trials';

            % helper function to filter cell arrays of strings, i.e. enable
            % selectors such as: filter_trials('start', {'east', 'west'})
            function mask = trial_filter(~, trial_data, selection)
                if ischar(selection)
                    mask = strcmp(trial_data,selection);
                elseif iscell(selection)
                    mask = false(size(trial_data));
                    for a = 1:length(selection)
                        mask = mask | strcmp(trial_data,selection{a});
                    end
                else
                    error('selection criterion must be a cell array or string')
                end
            end

        end   

        function highlight_cell(obj, cell_indices)
            obj.plot_cell_map({cell_indices, 'c'}, 'enable_class_colors');
        end
        
        % Accessors
        %------------------------------------------------------------
        function count = num_classified_cells(obj)
            count = sum(obj.is_cell);
        end
        
        function traces = get_trial(obj, trial_idx, varargin)
            % Return the traces for _all_ cells for the selected trial.
            %
            % TODO:
            %   - z-scoring option,
            %   - more precise baseline determination
            normalize_traces = false;
            
            for k = 1:length(varargin)
                vararg = varargin{k};
                switch lower(vararg)
                    case 'norm'
                        normalize_traces = true;
                end
            end
            
            traces = obj.trials(trial_idx).traces;
            if normalize_traces
                for j = 1:obj.num_cells
                    min_j = obj.trace_range(j,1);
                    max_j = obj.trace_range(j,2);
                    traces(j,:) = (traces(j,:)-min_j)/(max_j-min_j);
                end
            end
        end
        
        function [trace, frame_indices, selected_trials] = get_trace(obj, cell_idx, varargin)
            % For a cell, return traces over selected sets of trials.
            normalize_trace = false;
            selected_trials = 1:obj.num_trials;
            fill_type = 'traces';
            
            for k = 1:length(varargin)
                vararg = varargin{k};
                if isnumeric(vararg)
                    selected_trials = vararg;
                elseif ischar(vararg)
                    switch lower(vararg)
                        % TODO: Explore complications (if any) between
                        % normalization and trace fill type
                        case 'norm'
                            normalize_trace = true;
                        case 'fill'
                            fill_type = varargin{k+1};
                    end
                end
            end
            filtered_trials = find(obj.filter_trials(varargin{:}));
            selected_trials = intersect(selected_trials, filtered_trials)';
            
            trace = [];
            frame_indices = [];
            for k = selected_trials
                tr = obj.trials(k).traces(cell_idx,:);
                if ~obj.is_eventdata_loaded
                    % If eventdata is not available, then the only fill
                    % type (currently) allowed is 'traces'
                    trf = tr;
                else
                    ed = obj.trials(k).events{cell_idx};
                    trf = zeros(size(tr));
                    switch fill_type
                        case {'trace', 'traces'}
                            trf = tr;

                        case 'copy'                       
                            for m = 1:size(ed,1)
                                ef = ed(m,1):ed(m,2); % event frames (trough to peak)
                                trf(ef) = tr(ef);
                            end

                        case {'copyamp', 'copyzero'}
                            for m = 1:size(ed,1)
                                ef = ed(m,1):ed(m,2);
                                trf(ef) = tr(ef) - tr(ef(1));
                            end

                        otherwise
                            error('Fill type "%s" not recognized', fill_type);
                    end
                end
                trace = [trace trf]; %#ok<*AGROW>
                frame_indices = [frame_indices obj.trial_indices(k,1):obj.trial_indices(k,end)];
            end
            
            if normalize_trace
                trace_min = obj.trace_range(cell_idx,1);
                trace_max = obj.trace_range(cell_idx,2);
                trace = (trace - trace_min) / (trace_max - trace_min);
            end
        end
        
        function [traces, align_info] = get_aligned_trace(obj, cell_idx, trial_inds, alignment_frames, varargin)
            % Aligns a single cell trace across trials, and returns the
            % result as a [num_trials x num_common_frames] matrix.

            if ~isempty(varargin)
                for k = 1:length(varargin)
                    vararg = varargin{k};
                    if ischar(vararg)
                        switch lower(vararg)
                            % If the provided 'alignment_frames' is
                            % computed with respect to the start of each
                            % trial, then trial offsets need to be applied
                            % before use by 'get_aligned_trace'
                            case 'apply_trial_offset'
                                alignment_frames = alignment_frames +...
                                    (obj.trial_indices(trial_inds,1)-1);
                        end
                    end
                end
            end
            
            if islogical(trial_inds)
                trial_inds = find(trial_inds);
            end
            N = length(trial_inds);
            [pre_offset, post_offset] = compute_frame_offsets(obj.trial_indices, trial_inds, alignment_frames);
            num_common_frames = post_offset-pre_offset+1;
            
            traces = zeros(N, num_common_frames);
            
            for k = 1:N
                trial_ind = trial_inds(k);
                af = alignment_frames(k) - (obj.trial_indices(trial_ind,1)-1);
                pre_frame = af + pre_offset;
                post_frame = af + post_offset;
                tr = obj.get_trace(cell_idx, trial_ind, varargin{:});
                traces(k,:) = tr(pre_frame:post_frame);
            end
            
            align_info.num_trials = N;
            align_info.trial_inds = trial_inds;
            align_info.alignment_frames = alignment_frames;
            align_info.aligned_time = pre_offset:post_offset;
        end
        
        function es = get_events(obj, cell_idx, trial_idx, varargin)
            % Default options
            align_to = 1; % One of 1:4, indicating the frame within trial for alignment
            
            if ~isempty(varargin)
                for k = 1:length(varargin)
                    vararg = varargin{k};
                    switch lower(vararg)
                        case {'align', 'align_to'}
                            align_to = varargin{k+1};
                    end
                end
            end
            
            % Syntactic sugar for the following...
            es = obj.trials(trial_idx).events{cell_idx};
            
            if ~isempty(es)
                alignment_offset = diff(obj.trial_indices(trial_idx, [1 align_to]));
                es(:,1:2) = es(:,1:2) - alignment_offset;
            end
        end
        
        function es = get_events_full(obj, cell_idx)
            % Retrieves the set of events for all trials, adding in frame
            % offsets for each trial
            es = [];
            for k = 1:obj.num_trials
                es_k = obj.trials(k).events{cell_idx};
                if ~isempty(es_k)
                    trial_offset = obj.trial_indices(k,1) - 1;
                    es_k(:,1:2) = es_k(:,1:2) + trial_offset;
                    es = cat(1, es, es_k);
                end
            end
        end
        
        function mask = get_mask(obj, cell_indices)
            % When 'cell_indices' is omitted, then return the masks of all
            % classified cells
            if ~exist('cell_indices', 'var')
                cell_indices = find(obj.is_cell);
            end
            
            [height, width] = size(obj.cells(1).im);
            mask = zeros(height, width);
            for cell_idx = cell_indices
                mask = mask | obj.cells(cell_idx).mask;
            end
        end
              
        function is_correct = get_trial_correctness(obj)
            is_correct = [obj.trials.correct];
        end
        
        function selected_idx = get_cell_by_xy(obj, xy, varargin)
            % Returns the first cell index whose filter boundary encloses
            % the XY position provided in 'xy'
            
            % By default all cell candidates can be clicked
            classified_cells_only = 0;
            for i = 1:length(varargin)
                if ischar(varargin{i})
                    switch lower(varargin{i})
                        case {'cells', 'cellsonly'}
                            classified_cells_only = 1;
                    end
                end
            end
            
            selected_idx = []; % If no hit, then return empty
            for k = 1:obj.num_cells
                boundary = obj.cells(k).boundary;
                if (obj.is_cell(k) || ~classified_cells_only)
                    if inpolygon(xy(1), xy(2), boundary(:,1), boundary(:,2))
                        selected_idx = k;
                        break;
                    end
                end
            end
        end
        
        function neighbor_inds = get_nearest_sources(obj, cell_idx, num_neighbors, varargin)
            % Return the indices of 'num_neighbors' number of sources
            % closest to the source with 'cell_idx'.
            
            classified_cells_only = 0;
            for i = 1:length(varargin)
                if ischar(varargin{i})
                    switch lower(varargin{i})
                        case {'cell', 'cells', 'cellsonly'}
                            classified_cells_only = 1;
                    end
                end
            end

            d = obj.cell_distances(cell_idx,:); % Distance to all cells
            if classified_cells_only
                is_not_cell = ~obj.is_cell;
                d(is_not_cell) = Inf; % Set distances to non-cells as Inf
            end
            [~, neighbor_inds] = sort(d); % Ascending order
            neighbor_inds = neighbor_inds(1:num_neighbors);
        end
        
        % Classification
        %------------------------------------------------------------
        function class = get_class(obj)
            class = {obj.cells.label}'; % Column cell
        end
        
        function apply_labels_to(obj, label, cell_indices)
            for k = cell_indices
                obj.cells(k).label = label;
            end
        end
        
        function set_labels(obj, varargin)
            cell_indices = 1:obj.num_cells;
            if ~isempty(varargin)
                cell_indices = varargin{1};
            end
            obj.apply_labels_to('cell', cell_indices);
        end
        
        function reset_labels(obj, varargin)
            cell_indices = 1:obj.num_cells;
            if ~isempty(varargin)
                cell_indices = varargin{1};
            end
            obj.apply_labels_to([], cell_indices);
        end
        
        function set_unlabeled_cells(obj)
            unlabeled_cells = find(cellfun(@isempty, obj.get_class))';
            obj.apply_labels_to('not a cell', unlabeled_cells);
        end
        
        function invert_labels(obj)
            orig_cells = find(obj.is_cell);
            orig_not_cells = setdiff(1:obj.num_cells, orig_cells);
            obj.apply_labels_to('not a cell', orig_cells);
            obj.apply_labels_to('cell', orig_not_cells);
        end
        
        function is_cell = is_cell(obj, cell_indices)
            % When 'cell_indices' is omitted, then return the label of all
            % cells
            if ~exist('cell_indices', 'var')
                cell_indices = 1:obj.num_cells;
            end
            
            is_cell = zeros(size(cell_indices));
            for k = 1:length(cell_indices)
                cell_idx = cell_indices(k);
                is_cell(k) = any(strcmp(obj.cells(cell_idx).label,...
                    {'phase-sensitive cell', 'cell'}));
            end
        end
        
        % Load behavior movie
        %------------------------------------------------------------
        function loaded = is_behavior_loaded(obj)
            loaded = ~isempty(obj.behavior_vid);
        end
        
        function load_behavior_movie(obj, behavior_source)
            obj.behavior_vid = VideoReader(behavior_source);
            fprintf('%s: Loaded behavior video from "%s"\n',...
                datestr(now), behavior_source);
            if (obj.behavior_vid.NumberOfFrames ~= obj.full_num_frames)
                fprintf('  Warning! Number of frames in behavior video (%d) does not match the trial frame table (%d)!\n',...
                    obj.behavior_vid.NumberOfFrames, obj.trial_indices(end,end));
            end
            
            % Load reference image
            img = obj.behavior_vid.read(1);
            obj.behavior_ref_img = squeeze(img(:,:,1,:));
        end
        
        function Mb = get_behavior_trial(obj, trial_idx)
            % Must use original frame indices to access raw data
            trial_frames = obj.orig_trial_indices(trial_idx, [1 end]); % [Start end]
            Mb = obj.behavior_vid.read(trial_frames);
            Mb = squeeze(Mb(:,:,1,:)); % Movie is actually grayscale!
        end
        
        function A = get_behavior_trial_frame(obj, trial_idx, trial_frame_idx)
            frame_idx = obj.orig_trial_indices(trial_idx, 1) + trial_frame_idx - 1;
            A = obj.behavior_vid.read(frame_idx);
            A = squeeze(A(:,:,1));
        end
        
        % Load tracking data
        %------------------------------------------------------------
        function load_tracking(obj, tracking_source)
            centroids = load(tracking_source);

            % Must use original frame indices to access raw data
            for k = 1:obj.num_trials
                inds = obj.orig_trial_indices(k,1):obj.orig_trial_indices(k,end);
                obj.trials(k).centroids = centroids(inds, :);
            end
            
            fprintf('%s: Loaded tracking data from "%s"\n',...
                datestr(now), tracking_source);
            obj.is_tracking_loaded = true;
        end
        
        % Load event data
        %------------------------------------------------------------
        function load_events(obj, event_source)
            data = load(event_source);
            assert(length(data.events) == obj.num_cells,...
                'Error: Number of cells in event file does not match that in DaySummary!');
            
            % Note: We are expecting that event detection has been run on
            % _all_ trials, including probes.
            assert(data.events(1).info.num_frames == obj.full_num_frames,...
                'Error: Number of frames in event file does not match full number of frames in DaySummary!');
            
            events_per_trial = compute_events_per_trial({data.events.auto},...
                obj.orig_trial_indices, obj.full_num_frames);
            for k = 1:obj.num_trials
                obj.trials(k).events = events_per_trial{k};
            end
            
            fprintf('%s: Loaded events from "%s"\n',...
                datestr(now), event_source);
            obj.is_eventdata_loaded = true;
        end
        
        % Inspect switch data
        %------------------------------------------------------------
        function loaded = is_switchdata_loaded(obj)
            valid_starts = {'east', 'west'};
            loaded = ~isempty(obj.switchdata.pre_switch_trials) &&...
                     ~isempty(obj.switchdata.post_switch_trials) &&...
                     any(strcmp(obj.switchdata.changing_path_start, valid_starts)) &&...
                     any(strcmp(obj.switchdata.constant_path_start, valid_starts));
        end
        
        function st = get_switch_trials(obj)
            % Provide trial indices for switch analysis. Note that we only
            % examine correct trials.
            %
            % NOTE: Consider setting once when switchdata is loaded
            st = struct('constant_pre', [],...
                        'constant_post', [],...
                        'changing_pre', [],...
                        'changing_post', []);
                    
            if obj.is_switchdata_loaded
                sd = obj.switchdata;
                st.constant_pre = find(obj.filter_trials(...
                    'range', sd.pre_switch_trials,...
                    'start', sd.constant_path_start,...
                    'correct'));
                st.constant_post = find(obj.filter_trials(...
                    'range', sd.post_switch_trials,...
                    'start', sd.constant_path_start,...
                    'correct'));
                st.changing_pre = find(obj.filter_trials(...
                    'range', sd.pre_switch_trials,...
                    'start', sd.changing_path_start,...
                    'correct'));
                st.changing_post = find(obj.filter_trials(...
                    'range', sd.post_switch_trials,...
                    'start', sd.changing_path_start,...
                    'correct'));
            end
        end
        
        function trial_inds = get_constant_path_trials(obj)
            trial_inds = find(obj.filter_trials('start', obj.switchdata.constant_path_start));
        end
        
        function trial_inds = get_changing_path_trials(obj)
            trial_inds = find(obj.filter_trials('start', obj.switchdata.changing_path_start));
        end
        
        function reset_switchdata(obj)
            obj.switchdata.pre_switch_trials = [];
            obj.switchdata.post_switch_trials = [];
            obj.switchdata.changing_path_start = '';
            obj.switchdata.constant_path_start = '';
            obj.switchdata.changing_path_cutoff = [];
        end
            
    end % public methods

end
