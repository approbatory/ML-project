function X_all = gen_all_X_at_pos_closest(ds, poss)
    [frame_of_interest, ~] = find_frame_at_pos(ds, poss);

    X_all = zeros(ds.num_cells, ds.num_trials, size(poss, 2));
    for  j = 1:size(poss, 2)
        X_all(:,:,j) = gen_X_at_frames( {ds.trials.events}, frame_of_interest(:,j));
    end
end
