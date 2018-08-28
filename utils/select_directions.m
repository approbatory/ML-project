function [fw, bw] = select_directions(XY, varargin)
p = inputParser;
p.addRequired('XY', @(x) isnumeric(x) && size(x,2) == 2);
p.addParameter('old_method', false, @islogical);
p.addParameter('ranges', false, @islogical);
p.parse(XY, varargin{:});
XY = p.Results.XY;
old_method = p.Results.old_method;
ranges = p.Results.ranges;
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
    
    if trial_end(1) <= trial_start(1)
        trial_end = trial_end(2:end);
    end
    
    fw = false(length(track_coord), 1);
    bw = false(length(track_coord), 1);
    %%%
    fw_ranges.start = [];
    fw_ranges.end = [];
    bw_ranges.start = [];
    bw_ranges.end = [];
    %%%
    for tr_i = 1:min(length(trial_start),length(trial_end))
        trial_vel = velocity(trial_start(tr_i):trial_end(tr_i));
        
        if all(trial_vel > vthresh)
            fw(trial_start(tr_i):trial_end(tr_i)) = true;
            fw_ranges.start = [fw_ranges.start trial_start(tr_i)];
            fw_ranges.end = [fw_ranges.end trial_end(tr_i)];
        elseif all(trial_vel < -vthresh)
            bw(trial_start(tr_i):trial_end(tr_i)) = true;
            bw_ranges.start = [bw_ranges.start trial_start(tr_i)];
            bw_ranges.end = [bw_ranges.end trial_end(tr_i)];
        end
    end
    if ranges
        fw = fw_ranges;
        bw = bw_ranges;
    end
    
end
end