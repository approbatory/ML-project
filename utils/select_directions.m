function [fw, bw] = select_directions(XY)
if isstruct(XY)
    XY = XY.trials.centroids;
end
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
end