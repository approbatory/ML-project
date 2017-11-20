function X_all = gen_all_X_at_pos_closest_shuffled_start_keeped(ds, poss)

    X_all = gen_all_X_at_pos_closest(ds, poss);
    incorrect_trials = cell(1,4);
    correct_trials = cell(1,4);
    arms = ["east", "west", "south", "north"];
    
    for i =1:4
      incorrect_trials{i} = find(cell2mat({ds.trials.correct}) == false & ...
                                 string({ds.trials.start}) == arms(i) );
      correct_trials{i} = find(cell2mat({ds.trials.correct}) == true & ...
                               string({ds.trials.start}) == arms(i) );
    end
    
    for n_pos = 1:size(poss,2)
        for n_cell = 1:ds.num_cells
           for i =1:4
                X_all(n_cell, incorrect_trials{i}, n_pos) = ...
                    X_all(n_cell, incorrect_trials{i}(randperm(size(incorrect_trials{i}, 2))), n_pos);

                X_all(n_cell, correct_trials{i}, n_pos) = ...
                    X_all(n_cell, correct_trials{i}(randperm(size(correct_trials{i}, 2))), n_pos);
            end
        end
    end
end
