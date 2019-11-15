function out_filter = threshold_ic_filter(in_filter, thresh)

M = max(in_filter(:));

out_filter = in_filter;
out_filter(out_filter < thresh*M) = 0;
out_filter = out_filter / sum(out_filter(:));