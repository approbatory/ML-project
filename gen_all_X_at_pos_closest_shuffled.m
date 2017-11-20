function X_all = gen_all_X_at_pos_closest_shuffled(ds, poss)

    X_all = gen_all_X_at_pos_closest(ds, poss);
    
    incorrect_trials = find(cell2mat({ds.trials.correct}) == false);
    correct_trials = find(cell2mat({ds.trials.correct}) == true);
    
    for n_pos = 1:size(poss,2)
        for n_cell = 1:ds.num_cells
           
            X_all(n_cell, incorrect_trials, n_pos) = ...
                X_all(n_cell, incorrect_trials(randperm(size(incorrect_trials, 2))), n_pos);

            X_all(n_cell, correct_trials, n_pos) = ...
                X_all(n_cell, correct_trials(randperm(size(correct_trials, 2))), n_pos);
        end
    end
end
