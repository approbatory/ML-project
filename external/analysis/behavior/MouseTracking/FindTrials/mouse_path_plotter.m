function [trial_frames] = mouse_path_plotter( movie, centroids_xy, varargin)
% Plots the path of the mouse during each trial based on position data and
% trial frames (either specified by plusmaze txt or extracted from movie
% using find_start_end_of_trials). Colormap = Jet --> Blue-ish = start of
% trial, Red-ish = end of trial
% 
% Input:
%     movie: Behavior video
%     centroids_xy: position data (.xy file) from get_mouse_XY_pos
%     varargin: plusmaze txt file, e.g. 'mouse7_day09_allo-south.txt'
% 
% Output:
%     trial_frames: if no plusmaze source, can return trial_frames generated
%         by find_start_end_of_trials
% 
% 2015-03-01 Fori Wang
% 2016-04-11 Updated Fori Wang

    % read in plusmaze source if available
    if ~isempty(varargin)
        maze_data = varargin{1};
        [frame_indices,~,~]=parse_plusmaze(maze_data);
        trial_frames = [frame_indices(:,1) frame_indices(:,4)];
    else
        % if no plusmaze source, find trial start and ends from movie
        [trial_frames,~]=find_start_end_of_trials(centroids_xy);
    end
    
    % load centroids
    centroids = load(centroids_xy);
    
    % read in behavior movie
    fprintf('Reading in Behavior Movie...\n');
    behavior_vid = VideoReader(movie);
    fprintf('Done!\n');
    fprintf('Plotting mouse paths...\n');
    figure;
    
    
    for trial_idx = 1:size(trial_frames,1)
        trial_start = trial_frames(trial_idx,1);
        trial_end = trial_frames(trial_idx,2);
        scrollsubplot(5,6,trial_idx)
        
        %set image as first image of trial
        image = read(behavior_vid,trial_start);
        imagesc(image);
        title(sprintf('Trial %d',trial_idx));
        axis image;
        set(gca,'XTick',[],'YTick',[]);
        colormap gray;
        hold on
        
        color_options = jet(size(centroids(trial_start:trial_end,1),1));
        % plot centroids in color
        for i = trial_start:trial_end
            plot(centroids(i,1),centroids(i,2),'*','Color',color_options(i-trial_start+1,:))
        end
    end

end