function trf = filter_trace(tr, cutoff_freq, sampling_freq)
% Zero-phase filtering by 2nd-order Butterworth lowpass filter

[b,a] = butter(2, cutoff_freq/(sampling_freq/2));
trf = filtfilt(b,a,double(tr));