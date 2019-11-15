function [X, xs, ys] = export_traces_naive(md, trial_map, extent)
% X = EXPORT_TRACES_NAIVE(md, trial_map, extent)
%
% Exports traces of all trials specified by trial_map into a 3D matrix 'X'
%   with dimensions [neurons x time x trials]. Each trial has been time
%   normalized to have the same length (in number of samples). Parameter
%   'extent' is one of 'full', 'first', 'second' and allows for
%   specification of sub-portions of trials.
%
% Also provides the mouse's trajectory on each trial. Trajectories are also
%   time normalized.
%

    num_cells = md.num_cells;
    num_samples = compute_resample_size(md); % Num samples per trial
    num_trials = size(trial_map,1);
    
    X = zeros(num_cells, num_samples, num_trials);
    
    xs = cell(num_trials,1);
    ys = cell(num_trials,1);
    
    if ~strcmp(extent, 'full')
        assert(all_tracking_loaded(md),...
               'Use of extent other than "full" requires tracking data for each day of MultiDay!');
    end
    
    resample_grid = linspace(0, 1, num_samples);
    for k = 1:num_trials
        % day and neuron indices
        di = trial_map(k,1);
        ni = md.get_indices(di);

        % traces for this trial
        trial = md.day(di).trials(trial_map(k,2));
        [t,x,y] = truncate_trial(trial.centroids, trial.start, extent);
        traces = trial.traces(ni,t);
        
        % Resample each trial to common number of samples
        orig_grid = linspace(0, 1, size(traces,2));
        X(:,:,k) = interp1(orig_grid, traces', resample_grid)';
        
        if md.day(di).is_tracking_loaded
            xs{k} = interp1(orig_grid, x, resample_grid);
            ys{k} = interp1(orig_grid, y, resample_grid);
        end
    end

    function [t_idx,x,y] = truncate_trial(xy, start_arm, extent)
        % x and y coordinates
        x = xy(:,1);
        y = xy(:,2);

        % truncate as specified
        switch lower(extent)
            case 'full'
                t_idx = true(size(x));
            
            case 'first'
                switch start_arm
                    case 'east'
                        t_idx = y > east_start_boundary(x);
                    case 'west'
                        t_idx = y < west_start_boundary(x);
                    otherwise
                        error('extent not implemented for probe trials');
                end
                % For 'first' half of trial, we expect t_idx to be 1 at the
                % beginning of the trial, and 0 at end. Filter for the FIRST
                % contiguous block of 1s.
                t_idx = cumprod(t_idx);
            
            case 'second'
                switch start_arm
                    case 'east'
                        t_idx = y < east_start_boundary(x);
                    case 'west'
                        t_idx = y > west_start_boundary(x);
                    otherwise
                        error('extent not implemented for probe trials');
                end
                % For 'second' half of trial, we expect t_idx to be 0 at the
                % beginning of the trial, and 1 at the end. Filter for the LAST
                % contiguous block of 1s.
                t_idx = flipud(cumprod(flipud(t_idx)));
            
            otherwise
                error('extent not specified correctly')
        end
        
        t_idx = logical(t_idx);
        x = x(t_idx,:);
        y = y(t_idx,:);

    end % truncate_trial

    % FIXME: Hard-coded trajectory parameters!
    function y = east_start_boundary(x)
        y = -x+600;
    end

    function y = west_start_boundary(x)
        y = -x+450;
    end
end % export_traces_naive

function all_loaded = all_tracking_loaded(md)
    % Check that each day of MultiDay has associated tracking data
    is_tracking_loaded = zeros(1, md.num_days);
    for k = 1:md.num_days
        is_tracking_loaded(k) = md.day(md.valid_days(k)).is_tracking_loaded;
    end
    all_loaded = all(is_tracking_loaded);
end % all_tracking_loaded

function num_samples = compute_resample_size(md)
% Compute the common number of samples to use when resampling traces.
    trial_indices = [];
    for di = md.valid_days
        trial_indices = [trial_indices; md.day(di).trial_indices]; %#ok<AGROW>
    end
    num_frames_per_trial = trial_indices(:,4) - trial_indices(:,1) + 1;

    % THK: Is max the right number of samples? Does it matter?
    num_samples = max(num_frames_per_trial);
end % compute_resample_size