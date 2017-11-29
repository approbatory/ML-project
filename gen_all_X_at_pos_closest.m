function X_all = gen_all_X_at_pos_closest(ds, poss)
    [frame_of_interest, pos] = find_frame_at_pos(ds, poss);

    X_all = zeros(ds.num_trials, ds.num_cells, length(poss));
    for  j = 1:length(poss)
        X_all(:,:,j) = gen_X_at_frames( {ds.trials.events}, frame_of_interest(:,j));
    end
end
