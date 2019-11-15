function trace_out = fix_baseline(trace_in)
% Remove baseline offset of cell trace with following steps:
%    - sample baseline for overlapping small segments
%    - smooth baseline samples
%    - interpolate samples to 1:length(trace) to get the time varying
%      baseline estimate
%    - subtract baseline estimate from the trace

    k = 100; % Sampling period
    s = 1000; % Sample size
    c = 0.15; % scale for the smoothing filter length
    
    if size(trace_in,1)==1,trace_in=trace_in';end % make input a column vector
    
    % Sample baseline at uniformly spaced points
    num_frames = length(trace_in);
    num_samples = floor(num_frames/k);
    baselines = zeros(num_samples,2);
    for idx_sample = 0:num_samples-1    
        idx_begin = max(1,1+idx_sample*k-s/2);
        idx_end = min(num_frames,idx_sample*k+s/2);
        tr = trace_in(idx_begin:idx_end);
        % Baseline estimate is the most frequent bin in the histogram
        [counts,centers] = hist(tr,20);
        baseline = centers(find(counts==max(counts),1));
        baselines(idx_sample+1,:) = [1+idx_sample*k,baseline];
    end

    % Smooth the baseline values with moving average filter
    filt_hlen = floor(num_samples*c/2);% Half-length of the filter
    filt_out = zeros(num_samples,1);
    for idx_sample = 0:num_samples-1 % Manual convolution (due to edge effects)   
        idx_begin = max(1,1+idx_sample-filt_hlen);
        idx_end = min(num_samples,idx_sample+filt_hlen+1);
        b = baselines(idx_begin:idx_end,2);
        filt_out(idx_sample+1) = mean(b);
    end

    % Interpolate the sample baseline points
    subt = interp1(baselines(:,1),filt_out,(1:num_frames)','spline');

    trace_out = trace_in-subt;

 
