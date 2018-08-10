function [fw, bw] = select_directions(XY, old_method)
if nargin == 1
    old_method = false;
end
if isstruct(XY)
    XY = XY.trials.centroids;
end
if old_method
    %vthresh = 15 / 20; %~3cm/s
    vthresh = 4 .* max(max(XY)-min(XY)) ./ 120 ./ 20; %4 cm/s
    track_coord = XY(:,1);
    velocity = diff(track_coord);
    smooth_velocity = medfilt1(velocity, 21);
    is_forward = smooth_velocity > vthresh;
    is_backward = smooth_velocity < -vthresh;
    
    track_range = range(track_coord);
    track_min = min(track_coord);
    track_max = max(track_coord);
    bottom_twelfth = track_min + track_range/12;
    top_twelfth = track_max - track_range/12;
    
    in_between = (track_coord > bottom_twelfth) & (track_coord < top_twelfth);
    fw = [false; is_forward] & in_between;
    bw = [false; is_backward] & in_between;
else
    vthresh = 4 .* max(max(XY)-min(XY)) ./ 120 ./ 20; %4 cm/s
    track_coord = XY(:,1);
    velocity = [0; diff(track_coord)];
    %smooth_velocity = medfilt1(velocity, 21);
    
    track_range = range(track_coord);
    track_min = min(track_coord);
    track_max = max(track_coord);
    bottom_twelfth = track_min + track_range/12;
    top_twelfth = track_max - track_range/12;
    
    in_between = (track_coord > bottom_twelfth) & (track_coord < top_twelfth);
    trial_start = find(diff(in_between) == 1);
    trial_end = find(diff(in_between) == -1);
    
    
    fw = false(length(track_coord), 1);
    bw = false(length(track_coord), 1);
    for tr_i = 1:length(trial_start)
        trial_vel = velocity(trial_start(tr_i):trial_end(tr_i));
        
        if all(trial_vel > vthresh)
            fw(trial_start(tr_i):trial_end(tr_i)) = true;
        elseif all(trial_vel < -vthresh)
            bw(trial_start(tr_i):trial_end(tr_i)) = true;
        end
    end
    
    
end
end