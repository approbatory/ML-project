function [m,err] = agg_lists(list, groups)

assert(all(cellfun(@isvector, list)));

lengths = cellfun(@numel, list);
max_2 = maxk(lengths,2);
top = max_2(2);

[m,err] = deal(zeros(top,1));
for i = 1:top
    filt = lengths >= i;
    
    filt_list = list(filt);
    filt_groups = groups(filt);
    vals = cellfun(@(x)x(i), filt_list);
    [m(i), err(i)] = grouped_mean(vals(:), filt_groups(:));
end