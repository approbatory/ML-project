function [idx_mv_onset,idx_turn_onset] = calc_behavioral_event_indices(cents,extrema_x,extrema_y,start_str,end_str)
% Calculates the frames of movement and turn onsets given the centroids of a trial.
%   extrema_x (extrema_y) contains the minimum and maximum x (y) values of the centroids in
%   a session.
%   start_str and end_str are the start and end arms.
%   if a particular index is not found, then 0 is returned in its place.
%

    % Constants
    angle_thresh = 25;
    pix_veloc_thresh = 2;
    
    num_frames = size(cents,1);
    
    % Get endpoint estimates
    west = [extrema_x(1),extrema_y(1)];
    east = [extrema_x(2),extrema_y(2)];
    south = [extrema_x(1),extrema_y(2)];
    north = [extrema_x(2),extrema_y(1)];

    % Find the two principal directions (going from east and south)
    ew_dir = west-east;
    ew_dir = ew_dir/norm(ew_dir);
    sn_dir = north-south;
    sn_dir = sn_dir/norm(sn_dir);

    % Estimate the movement/turn onset by displacement along the principal
    % directions
    switch(start_str)
        case 'east'
            start_dir = ew_dir;
        case 'west'
            start_dir = -ew_dir;
        case 'south'
            start_dir = sn_dir;
        case 'north'
            start_dir = -sn_dir;
    end

    switch(end_str)
        case 'east'
            end_dir = -ew_dir;
        case 'west'
            end_dir = ew_dir;
        case 'south'
            end_dir = -sn_dir;
        case 'north'
            end_dir = sn_dir;
    end

    disp_st = smooth(cents*start_dir',3);
    velocity_st = diff(disp_st);
    disp_end = smooth(cents*end_dir',3);  
    velocity_end = diff(disp_end);
    angle = atan(velocity_end./velocity_st)*180/pi;

    % Find movement onset
    for idx_mv_onset = 51:(num_frames-51)
        % Need 9 consecutive positive displacements
        if prod(velocity_st(idx_mv_onset:idx_mv_onset+8)>pix_veloc_thresh)>0 
            break;
        end
    end    
    % No movement onset is found
    if idx_mv_onset== num_frames-51
        idx_mv_onset=0;
    end

    % Find turn onset
    for idx_turn_onset = max((idx_mv_onset+1),51):(num_frames-51)
        % Need 9 consecutive positive displacements + minimum displacement angle
        if prod(velocity_end(idx_turn_onset:idx_turn_onset+8)>pix_veloc_thresh)>0 ...
                && angle(idx_turn_onset)>angle_thresh 
            break;
        end
    end    
    % No turn onset is found
    if idx_turn_onset== num_frames-51
        idx_turn_onset=0;
    end

end
