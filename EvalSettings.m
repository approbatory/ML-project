classdef EvalSettings
    properties
        %boolean, provide adjacent bins or not?
        adjacent_bins = false
        
        %integer, can be 1,2,3. How many equal chunks of trials to provide.
        %This many functions need to be provided as well
        trial_split = 1
        %integer, how many random trial splits to use. Only used if
        %trial_split ~= 1
        trial_split_reps = 1
        
        
        %string, can be 'tn' or 'tbn'. Include all bins ('tbn') or map over
        %bins ('tn')
        %Must be 'tn' if adjacent_bins is true
        input_mode = 'tn'
        
        %Can be 'max' or {'value', integer} or {'range', integer}
        % 'max' uses all the neurons
        % {'value', integer} evaluates at the specified number of neurons,
        % if possible
        % {'range', integer} evaluates at a range of values starting at 1,
        % integer, 2*integer, 3*integer,... n*integer, max
        num_neurons = 'max'
        
        
        %Can be 'max' or {'value', integer} or {'range', integer}
        % 'max' uses all the trials
        % {'value', integer} evaluates at the specified number of trials,
        % if possible
        % {'range', integer} evaluates at a range of values starting at 1,
        % integer, 2*integer, 3*integer,... n*integer, max
        num_trials = 'max'
        
        %A (jagged) list of lists of parameters to pass to the functions.
        %The evaluation will loop over all parameter combinations
        params = {}
    end
    
    methods
        function o = EvalSettings(varargin)
            p = inputParser;
            p.addParameter('badj', false, @islogical);
            p.addParameter('tsplit', 1, @(x)isequal(x,1) || isequal(x,2) || isequal(x,3));
            p.addParameter('reps', 1, @isnumeric);
            p.addParameter('mode', 'tn', @(x)isequal(x,'tn')||isequal(x,'tbn'));
            p.addParameter('N', 'max', @EvalSettings.validate_num);
            p.addParameter('T', 'max', @EvalSettings.validate_num);
            p.addParameter('params', {}, @EvalSettings.validate_params);
            
            p.parse(varargin{:});
            r = p.Results;
            o.adjacent_bins = r.badj;
            o.trial_split = r.tsplit;
            o.trial_split_reps = r.reps;
            o.input_mode = r.mode;
            o.num_neurons = r.N;
            o.num_trials = r.T;
            o.params = r.params;
            
            %check for inconsistencies
            %if o.trial_split_reps > 1 && o.trial_split == 1
            %    error('Cannot use reps if there is no trial split');
            %end
            if o.adjacent_bins && strcmp(o.input_mode, 'tbn')
                error('Cannot use adjacent_bins and tbn mode at the same time');
            end
        end
        
        function l = num_neurons_list(o, dt)
            if ischar(o.num_neurons) && strcmp(o.num_neurons, 'max')
                l = dt.total_neurons;
            elseif iscell(o.num_neurons) && strcmp(o.num_neurons{1}, 'value')
                l = o.num_neurons{2};
                assert(l <= dt.total_neurons);
            elseif iscell(o.num_neurons) && strcmp(o.num_neurons{2}, 'range')
                dn = o.num_neurons{2};
                l = dn:dn:dt.total_neurons;
                if l(1) ~= 1
                    l = [1 l];
                end
                if l(end) ~= dt.total_neurons
                    l(end+1) = dt.total_neurons;
                end 
            else
                error('should not reach here');
            end
        end
        
        function l = num_trials_list(o, dt)
            if ischar(o.num_trials) && strcmp(o.num_trials, 'max')
                l = dt.n_min_trials;
            elseif iscell(o.num_trials) && strcmp(o.num_trials{1}, 'value')
                l = o.num_trials{2};
                assert(l <= dt.n_min_trials);
            elseif iscell(o.num_trials) && strcmp(o.num_trials{2}, 'range')
                dn = o.num_trials{2};
                l = dn:dn:dt.n_min_trials;
                %if l(1) ~= 1
                %    l = [1 l];
                %end
                if l(end) ~= dt.n_min_trials
                    l(end+1) = dt.n_min_trials;
                end 
            else
                error('should not reach here');
            end
        end
    end
    
    methods(Static)
        function res = validate_num(n)
            res = false;
            if ischar(n) && strcmp(n, 'max')
                res = true;
                return;
            end
            if iscell(n) && numel(n) == 2
                label = n{1};
                amount = n{2};
                
                if ~( strcmp(label, 'value') || strcmp(label, 'range') )
                    res = false;
                    return;
                end
                
                if isnumeric(amount) && isreal(amount) && isscalar(amount)...
                        && round(amount)==amount && amount > 0
                    res = true;
                    return;
                end
            end
        end
        
        function res = validate_params(p)
            res = false;
            if iscell(p)
                ok = cellfun(@isvector, p);
                if all(ok)
                    res = true;
                    return;
                end
            end
        end
    end
end