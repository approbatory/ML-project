function X_all = gen_all_X_at_pos_closest_shuffled(ds, poss)

    X_all = gen_all_X_at_pos_closest(ds, poss);

    for n_pos = 1:size(poss,2)
        for n_cell = 1:ds.num_cells
        X_all(n_cell, :, n_pos) = ...
            X_all(n_cell, randperm(ds.num_trials), n_pos);
        end
    end
end
