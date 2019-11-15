function [X, trial_meta, export_info] = export_tensor(md, varargin)
% [X, trial_meta, export_info] = EXPORT_TENSOR(md)
%
% Exports cross-day aligned cell traces from MultiDay into a lightweight
% format for further analysis (e.g. tensor analysis, decoding, etc.).
%
% Note that, export_tensor does NOT perform any normalization across days!
%
% Output format:
%
%   X: Cross-day traces formatted into a [neurons x time x trials] matrix.
%       Note that all trials will be formatted to the same length, using
%       one of several "timewarp" methods.
%
%   trial_meta: Metadata associated with each trial (e.g. start loc, etc.)
%
%   export_info: Struct containing additional info such as:
%       info.neuron_map: Matrix [neurons x num_days] that maps each 
%                        matched neuron to its per-day neuron index.
%
%       info.trial_map: Matrix [trials x 2] where trial_map(i,1) is the day 
%                       index and trial_map(i,2) is the trial index in the
%                       original day.
%

    timewarp_method = 'align_to';

    extent = 'full'; % Used with 'naive' time warping
    align_idx = 3; % Used with 'align' time warping
    
    for k = 1:length(varargin)
        vararg = varargin{k};
        if ischar(vararg)
            switch lower(vararg)
                case 'method'
                    timewarp_method = lower(varargin{k+1});
                case 'align_to' % Used with 'align' timewarp
                    align_to = varargin{k+1};
                    align_events = {'start', 'open', 'close', 'end'};
                    align_idx = find(strcmp(align_events, align_to));
                    assert(~isempty(align_idx),...
                           'Align_to must be one of: %s', strjoin(align_events, ', '));
                case 'extent' % 'naive' timewarping
                    extent = varargin{k+1}; % 'full', 'first', or 'second'
            end
        end
    end
       
    % get trial map (filtering out those specified)
    trial_map = md.filter_trials(varargin{:});

    export_info.timewarp_method = timewarp_method;
    export_info.trial_map = trial_map;
    export_info.neuron_map = md.matched_indices;
    
    % activity traces for each trial
    switch timewarp_method
        case 'naive'
            % 'naive' method will take each trial, and rescale time so
            % that each trial has the same number of samples.
            fprintf('  Exporting MD with NAIVE time warping method (extent="%s")...\n', extent);
            [X,x,y] = export_traces_naive(md, trial_map, extent);
            export_info.naive.extent = extent;
            
        case {'align', 'align_to'}
            % 'align' method will align each trial to one of four
            % intra-trial events:
            %   (1) Start of trial
            %   (2) Opening of gate
            %   (3) Closing of gate
            %   (4) End of trial
            % taking as many samples that is common to all trials (i.e.
            % bounded by the shortest trial). Time is not scaled, and can
            % be directly compared between trials.
            fprintf('  Exporting MD with ALIGN time warping method (align_idx=%d)...\n', align_idx);
            [X, x, y, align_axis] = export_traces_align(md, trial_map, align_idx);
            export_info.align.idx = align_idx;
            export_info.align.axis = align_axis;

        otherwise
            error('Time warping method ("%s") not recognized.', timewarp_method);
    end
    
    % metadata for each trial (start, end, turn, correct, etc.)
    trial_meta = export_metadata(md, trial_map);
    trial_meta.x = x; % also export position
    trial_meta.y = y;
end % export
