function [ trial_frames, tracking_error_frames] = find_start_end_of_trials( centroids_xy )
% Using position data from get_mouse_XY_pos (using behavior vid), find
% frames where trial ends and next trial starts (based on fact that mouse
% will jump from e.g., south arm to east arm)
% 
%   ** Note: These trial boundaries must be confirmed by the user; in
%   egocentric sessions, probe trials will often be merged with trials
%   before or after (because there is no visible jump in position)
% 
% Input:
%     .xy file of position data from get_mouse_XY_pos.m
% 
% Returns
%     trial_frames: Array of start and end frames [start end; start end; ...]
%     tracking_error_frames: Frames that are not start/end but show a jump
% 
% 2015-02-28 Fori Wang
% Updated 2016-04-11 Fori Wang

    centroids = load(centroids_xy);
    num_frames = length(centroids);
    trial_counter = 1;
    tracking_error_frames = [];
    trial_frames = 1;
    
    for frame_idx = 2:num_frames
        
        % grab the current and previous centroid
        this_centroid = centroids(frame_idx,:);
        previous_centroid = centroids(frame_idx-1,:);
        x = this_centroid(1);
        y = this_centroid(2);
        x_old = previous_centroid(1);
        y_old = previous_centroid(2);
        
        % find jumps in position (potential start/end points)
        distance_thresh = 50;
        distance = sqrt((x-x_old)^2+(y-y_old)^2);
        
        if distance > distance_thresh
        
            % limit potential start/ends of trial to frames where start/end
            % is in one of the end arms
            west_x_thresh = 250;
            west_y_thresh = 150;
            east_x_thresh = 350;
            east_y_thresh = 300;
            south_x_thresh = west_x_thresh;
            south_y_thresh = east_y_thresh;
            north_x_thresh = east_x_thresh;
            north_y_thresh = west_y_thresh;

            % WEST start
            if x < west_x_thresh && y < west_y_thresh
                if x_old < south_x_thresh && y_old > south_y_thresh || ...% SOUTH end
                        x_old > north_x_thresh && y_old < north_y_thresh ||...% NORTH end
                        x_old < west_x_thresh && y_old < west_y_thresh ||...% WEST end
                        x_old > east_x_thresh && y_old > east_y_thresh % EAST end
                    trial_frames(trial_counter,2)=frame_idx-1; % end
                    trial_frames(trial_counter+1,1)=frame_idx; % next start
%                     trial_frames(trial_counter+1,3)=4; % mark start arm
                    trial_counter = trial_counter+1;
                end
            
            % EAST start
            elseif x > east_x_thresh && y > east_y_thresh
                if x_old < south_x_thresh && y_old > south_y_thresh || ...% SOUTH end
                        x_old > north_x_thresh && y_old < north_y_thresh ||...% NORTH end
                        x_old < west_x_thresh && y_old < west_y_thresh ||...% WEST end
                        x_old > east_x_thresh && y_old > east_y_thresh % EAST end
                    trial_frames(trial_counter,2)=frame_idx-1; % end
                    trial_frames(trial_counter+1,1)=frame_idx; % next start
%                     trial_frames(trial_counter+1,3)=2; % mark start arm
                    trial_counter = trial_counter+1;
                end
                
            % SOUTH start
            elseif x < south_x_thresh && y > south_y_thresh
                if x_old < south_x_thresh && y_old > south_y_thresh || ...% SOUTH end
                        x_old > north_x_thresh && y_old < north_y_thresh ||...% NORTH end
                        x_old < west_x_thresh && y_old < west_y_thresh ||...% WEST end
                        x_old > east_x_thresh && y_old > east_y_thresh % EAST end
                    trial_frames(trial_counter,2)=frame_idx-1; % end
                    trial_frames(trial_counter+1,1)=frame_idx; % next start
%                     trial_frames(trial_counter+1,3)=3; % mark start arm
                    trial_counter = trial_counter+1;
                end
                
            % NORTH start
            elseif x > north_x_thresh && y < north_y_thresh
                if x_old < south_x_thresh && y_old > south_y_thresh || ...% SOUTH end
                        x_old > north_x_thresh && y_old < north_y_thresh ||...% NORTH end
                        x_old < west_x_thresh && y_old < west_y_thresh ||...% WEST end
                        x_old > east_x_thresh && y_old > east_y_thresh % EAST end
                    trial_frames(trial_counter,2)=frame_idx-1; % end
                    trial_frames(trial_counter+1,1)=frame_idx; % next start
%                     trial_frames(trial_counter+1,3)=1; % mark start arm
                    trial_counter = trial_counter+1;                
                end

            else % Random jump from bad tracking
                tracking_error_frames = [tracking_error_frames; frame_idx-1 frame_idx];
            end
        end
    end
    trial_frames(trial_counter,2)=num_frames; % last frame = end of last trial
end

