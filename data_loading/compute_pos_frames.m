function pos_frames = compute_pos_frames(ds, poss)
    %finds the frame closest to specified poss based on the angle of the turn,
    %only within gate open-close times, including open time, excluding close
    %time, i.e. start 1, open 3, close 7, end 10, only use 3 4 5 6
    pos = preprocess_xy(ds);
    pos_frames = zeros(length(ds.trials), length(poss));
    reparam_pos = cellfun(@reparam, pos, 'UniformOutput', false);
    for i = 1:length(ds.trials)
        c = ds.trial_indices(i,:); s = c(2) - c(1) + 1; e = c(3) - c(1);
        r = reparam_pos{i}(s:e,:);
        [~, ix] = min((r - poss).^2);
        pos_frames(i,:) = int32(ix) + s - 1;
    end

    function r = reparam(r)
        r = 0.5-abs(0.5-r);
        r = atan2(r(:,2), r(:,1))/(pi/2);
        if r(end) < r(1)
            r = 1 - r;
        end
    end
end