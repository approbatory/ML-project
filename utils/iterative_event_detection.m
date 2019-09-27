function events = iterative_event_detection(traces, varargin)
p = inputParser;
p.addRequired('traces', @(x) size(x,1)>size(x,2));
p.addParameter('z_val', 3, @isscalar);
p.addParameter('fps', 20, @isscalar);
p.addParameter('kern_size', 8, @isscalar);
p.addParameter('med_window_size', 201, @isscalar);
p.parse(traces, varargin{:});

z_val = p.Results.z_val;

fps = p.Results.fps;

traces = traces - medfilt1(traces, p.Results.med_window_size);

events = zeros(size(traces));


for c_i = 1:size(traces,2)
    dff_trace = traces(:,c_i);
    
    cutoff_freq = 4/30 * fps;
    [b,a] = butter(2, cutoff_freq/(fps/2));
    
    tr_smooth = filtfilt(b,a,double(dff_trace));
    
    %figure; hold on;
    %plot(dff_trace);
    %plot(tr_smooth);
    
    
    peaks = false(size(tr_smooth));
    for i = 1:20
        %mu = mean(tr_smooth(~peaks));
        mu = prctile(tr_smooth(~peaks), 25);
        sigma = std(tr_smooth(~peaks));
        peaks_new = tr_smooth > mu + z_val*sigma;
        
        if isequal(peaks_new, peaks)
            peaks = peaks_new;
            %fprintf('broke at iteration %d\n', i);
            break
        end
        
        peaks = peaks_new;
    end
    
    peak_starts = [diff(peaks)==1;0];
    extended = conv(peak_starts, ones(p.Results.kern_size,1));
    extended(extended > 1) = 1;
    events(:,c_i) = extended(1:numel(peaks));
end
%%

%figure; hold on;
%plot(dff_trace);
%plot((tr_smooth - mu)./sigma);
%plot(peaks.*z_val);
%end