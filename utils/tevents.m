function [single_time, events_auto] = tevents(traces)
[num_cells, num_frames] = size(traces);
single_time = zeros(num_frames, num_cells);
fps = 20;
% We'll for events in a smoothed version of the trace
% Default parameters comes from cerebellar processing, where we used
%   - 30 Hz sampling frequency
%   - 4 Hz cutoff frequency
cutoff_freq = 4/30 * fps;
events_auto = cell(1,num_cells);
for cell_idx = 1:num_cells
    trace_orig = traces(cell_idx,:);
    %using func filter_trace
    trace = filter_trace(trace_orig, cutoff_freq, fps);
    %using func estimate_baseline_sigma
    [baseline, sigma, ~] = estimate_baseline_sigma(trace);
    info = struct('num_frames', num_frames,...
        'baseline', baseline,...
        'sigma', sigma,...
        'threshold', baseline + 5*sigma,...
        'amp_threshold', 0.1);
    %using func find_events_in_trials
    events_auto{cell_idx} = find_events_in_trials(trace, [1 2 3 num_frames],...
        info.threshold, info.baseline, info.amp_threshold);
    
    if ~isempty(events_auto{cell_idx})
        start_ = events_auto{cell_idx}(:,1);
        end_ = events_auto{cell_idx}(:,2);
        %midp{cell_idx} = floor(start_ + end_);
        %jays{cell_idx} = cell_idx + zeros(length(start_),1);
        single_time(floor((start_ + end_)/2),cell_idx) = 1;
    end
end