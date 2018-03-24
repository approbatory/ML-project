function X = gen_place_decoding_X(ds, copy_zeroed)
if nargin == 1
    copy_zeroed = false;
end
X = cell(ds.num_trials,1);
for i = 1:ds.num_trials
    X{i} = zeros(size(ds.trials(i).centroids,1),ds.num_cells);
    for j = 1:ds.num_cells
        for e = ds.trials(i).events{j}'
            if isinf(e(1))
                s = 1;
            else
                s = e(1);
            end
            if ~copy_zeroed
                X{i}(s:e(2),j) = e(3);
            else
                X{i}(s:e(2),j) = ds.trials(i).traces(j,s:e(2));
                if ~isinf(e(1))
                     X{i}(s:e(2),j) = X{i}(s:e(2),j) - ds.trials(i).traces(j,s);
                end
            end
        end
    end
end
end