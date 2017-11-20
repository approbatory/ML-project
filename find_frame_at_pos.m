function [frame_of_interest, pos] = find_frame_at_pos(ds, poss)
%FIND_FRAME_AT_POS For each position in poss, find frame closest to it in
%each trial
%   Return the relavant frames in each trial for each value in poss 'frame_of_interest'
%   and the preprocessed xy positions in 'pos'
pos = preprocess_xy(ds);
start_labels = {ds.trials.start};
frame_of_interest = zeros(length(ds.trials), length(poss));
for i = 1:length(ds.trials)
%     distance_to_center = sqrt(sum( (pos{i}-[0.5, 0.5]).^2, 2));
%     [M, I] = min(distance_to_center);
%     path_distace = [M - distance_to_center(1:(I-1)); distance_to_center(I:end) - M ];
%     tr_pos = (path_distace - min(path_distace))/(max(path_distace) - min(path_distace));

    xs = pos{i}(:,1);
    ys = pos{i}(:,2);
    sl = start_labels{i};
    if strcmp(sl, 'east')
        tr_pos = 1 - xs;
    elseif strcmp(sl, 'west')
        tr_pos = xs;
    elseif strcmp(sl, 'north')
        tr_pos = 1 - ys;
    elseif strcmp(sl, 'south')
        tr_pos = ys;
    else
        error('%s is an invalid value',sl);
    end
    start_to_end_indexes = ...
        (ds.trial_indices(i,2):ds.trial_indices(i,3)) - ds.trial_indices(i,1) + 1;
    
    for j = 1:length(poss)
        frame_after = find(tr_pos(start_to_end_indexes) > poss(j),1) ...
            + start_to_end_indexes(1) - 1;
        if isempty(frame_after)
            continue;
        end
        if (frame_after > 1) && (abs(tr_pos(frame_after-1)-poss(j)) < abs(tr_pos(frame_after)-poss(j)))
            frame_of_interest(i,j) = frame_after - 1;
        else
            frame_of_interest(i,j) = frame_after;
        end
    end
end
end