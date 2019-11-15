function idx_event_frames = extract_candidate_events(trace,varargin)
% Extract maximum points in the trace as candidate events
% varargin : 1 optional argument. Set it to 1 if you want to plot the
% output

num_frames = length(trace);
mad_scale = 9;
h = SGsmoothingfilter(20,6);
smooth_trace = conv(trace,h,'same');
mad = compute_mad(smooth_trace);
thresh = mad_scale * mad;
idx_event_frames = calc_localmax_above_threshold(thresh,smooth_trace);

%go back to the original trace and spot the maximum
for k = 1:length(idx_event_frames)
    coarse_event = idx_event_frames(k);
    go_back = min(10,coarse_event-1);
    go_forward = min(5,num_frames-coarse_event);
    [~,idx_max] = max(trace((coarse_event-go_back):(coarse_event+go_forward)));
    idx_event_frames(k) = idx_max+coarse_event-go_back-1;
end



if ~isempty(varargin)
    if varargin{1} == 1
        spikes_vec = zeros(1,length(smooth_trace));
        spikes_vec(idx_event_frames) = 1;
        plot(smooth_trace);
        hold on;
        stem(spikes_vec*thresh);
        plot(ones(1,length(smooth_trace))*thresh,'r')
        plot(trace,'--m');
        hold off
        for k = 1:length(idx_event_frames)
            text(idx_event_frames(k),double(thresh*1.1),num2str(k))
        end
    end
end