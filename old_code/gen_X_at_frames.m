function X = gen_X_at_frames(evs, frames)
%generates the feature vector from the given events at the given frames
%training examples are in rows
assert(length(evs)==length(frames));
n_trials = length(evs);
n_cells = length(evs{1});
X = zeros(n_trials, n_cells);
for i = 1:n_trials
    for j = 1:n_cells
        for e = evs{i}{j}'
            if frames(i) >= e(1) && frames(i) <= e(2)
                X(i,j) = e(3);
            end
        end
    end
end
end