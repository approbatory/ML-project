function events = find_events(trace, threshold, baseline)
% Find segments of 'trace' above the specified 'threshold'. Every local
% maximum within those segments are identified as events. This approach to
% event detection tends to work well with smoothed (e.g. low-pass filtered)
% versions of calcium traces.
%
% Returns:
%   events: [num_events x 3] where
%     events(k,1): Frame of the trough preceding the k-th event.
%     events(k,2): Frame corresponding to the peak of the k-th event
%     events(k,3): Amplitude difference between peak and trough
%

above_thresh_frames = find(trace >= threshold);

if ~any(above_thresh_frames)
    events = [];
    return;
end

segments = frame_list_to_segments(above_thresh_frames);
num_segments = size(segments,1);

eventpeaks = cell(1,num_segments);
for k = 1:num_segments
    seg = segments(k,1):segments(k,2);
    tr_seg = trace(seg);
    
    eventpeaks{k} = find_all_localmax(seg, tr_seg);
end
eventpeaks = cell2mat(eventpeaks);

% Format: [Trough-preceding-event Event-peak Event-amplitude]
num_events = length(eventpeaks);
events = zeros(num_events, 3);
for k = 1:num_events
    peak_frame = eventpeaks(k);
    
    % Attempt to find the trough preceding the event peak. Note that there
    % are cases where the trough cannot be found:
    %   1) The peak is at the beginning of the trace
    %   2) The trough is at the beginning of the trace -- which means that
    %   it's possible that we didn't roll all the way down to the trough
    trough_not_found = false;
    if peak_frame == 1
        trough_not_found = true;
    else
        trough_frame = seek_localmin(trace,peak_frame-1);
        if (trough_frame == 1)
            trough_not_found = true;
        end
    end
    
    % In the case that trough was not found, use the difference between the
    % baseline and the event peak as the amplitude
    if trough_not_found
%         trough_frame = -Inf; % Flag for indicating that trough wasn't found
        trough_frame = 1;
        event_amp = trace(peak_frame) - baseline;
    else
        event_amp = trace(peak_frame) - trace(trough_frame); 
    end
    
    events(k,:) = [trough_frame peak_frame event_amp];
end

end % find_events

function events = find_all_localmax(seg, tr_seg)
    % Whether the first and last points of the segment can be
    % identified as a local maximum
    allow_ends = true;
    x1 = [allow_ends tr_seg(2:end)>tr_seg(1:end-1)];
    x2 = [tr_seg(1:end-1)>tr_seg(2:end) allow_ends];

    events = seg(x1 & x2);
end % find_all_localmax