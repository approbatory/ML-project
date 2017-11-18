function [pos] = preprocess_xy(ds)
%PREPROCESS_XY Summary of this function goes here
%   Detailed explanation goes here
n_trials = length(ds.trials);

pos = cell(n_trials,1);
for i = 1:n_trials
    pos{i} = ds.trials(i).centroids * [1 1;1 -1]/sqrt(2);
end
all_pos = cell2mat(pos);

mins = min(all_pos);
maxs = max(all_pos);

%all_pos_resc = (all_pos - mins)./(maxs - mins);
for i = 1:n_trials
    pos{i} = (pos{i} - mins)./(maxs-mins);
end

end

