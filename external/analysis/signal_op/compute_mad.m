function mad = compute_mad(trace)

med = median(trace);
mad = median(abs(trace-med));