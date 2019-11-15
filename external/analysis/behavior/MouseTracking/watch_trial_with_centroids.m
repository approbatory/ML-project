function watch_trial_with_centroids( movie, centroids_xy, trial, varargin )
% Plays behavior video with centroid on top for specified trial
% 
% Input:
%     movie: behavior video
%     centroids_xy: position data from get_mouse_XY_pos
%     trial: matrix of trial(s) to watch, e.g. [1 4 5 8]
%     varargin: plusmaze txt file, e.g. 'mouse7_day09_allo-south.txt' or
%               trimmed/binned version
%
% Example use:
% trials = 1:15;
% watch_trial_with_centroids('mouse7.mp4','centroids.mat',trials,'mouse7.txt')
%
% 2015-03-06 Fori Wang
%
%

    % get trial indices
    if ~isempty(varargin)
        trial_indices = varargin{1};
        trial_indices = get_trial_frame_indices(trial_indices);
        trial_indices = [trial_indices(:,1) trial_indices(:,4)];
    else
        [trial_indices,~] = find_start_end_of_trials(centroids_xy);
    end

    % read in behavior movie
    fprintf('Reading in Behavior Movie...\n');
    behavior_vid = VideoReader(movie);
    
    % load centroids
    centroids = load(centroids_xy);
    
    % initialize plot
    figure;
    
    this_trial_indices = trial_indices(1,:);
    trial_start = this_trial_indices(1);  
    first_image = read(behavior_vid,trial_start);
    k = imagesc(first_image);
    title(sprintf('Trial %d',trial(1)),'FontSize',18);
    axis image; colormap gray;
    hold on
    l = plot(centroids(trial_start,1),centroids(trial_start,2),'g*');
    
    blank_image = zeros(size(first_image));
    
    for i = 1:length(trial)
        
        % update trial indices
        this_trial = trial(i);
        title(sprintf('Trial %d',this_trial));
        this_trial_indices = trial_indices(this_trial,:);
        trial_start = this_trial_indices(1);  
        
        % load trial video
        trial_video = read(behavior_vid,this_trial_indices);
        trial_video = squeeze(trial_video(:,:,1,:)); % 3D movie stack
       
        % plot remaining frames
        for frame_idx = 1:size(trial_video,3);
            image = trial_video(:,:,frame_idx);
            this_centroid = centroids(frame_idx+trial_start-1,:);

            set(k,'CData',image);
            
            if frame_idx == 1
                centroid_color = 'm';
                pause_time = 1.0; % inspect first frame for bleedthrough
            else
                centroid_color = 'g';
                pause_time = 0.1; % go through rest of movie fast
            end
            
            set(l,'XData',this_centroid(1),'YData',this_centroid(2),'color',centroid_color);
            pause(pause_time); % play through rest of movie faster
            
        end
        
        % display black screen at end of each trial
        set(k,'CData',blank_image);
        set(l,'color','k');
        pause(0.25);
    end

end
