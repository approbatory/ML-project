function [ centroids ] = get_mouse_XY_pos( behavior_source, varargin )
% Extract the X,Y positions ( centroids ) of the mouse from the behavior
%   video ( movie ), with live monitor option (specify 'displayTracking')
%  
% Example uses:
%   Normal mode:
%       centroids = get_mouse_XY_pos('c9m7d08_ti2-sub.mp4');
%   Live monitor mode:	
%       centroids = get_mouse_XY_pos('c9m7d08_ti2-sub.mp4','displayTracking');
%   Trial-aware mode:
%       centroids = get_mouse_XY_pos('c9m7d08_ti2-sub.mp4','trialAware','c9m7d08_cr_ti2.txt');
% 
% Input:
%   - behavior_source: Behavior video (e.g. mp4, m4v)
%   - varargin:
%     'displayTracking': Displays side-by-side of original movie and
%           processing movie, along with centroid on each. blue * means good
%           centroid value, red * means could not find mouse, using previous
%           centroid value, magenta * means temporary centroid assigned at
%           trial boundary (to be reassigned same as next centroid)
%     'trialAware','trial_indices.txt': If 'trialAware' and .txt file provided,
%           will not use prev_centroid method at the trial boundaries to avoid
%           bleed-through across trials (i.e. taking the position of the
%           last frame -- in a previous trial -- when mouse is not detected).
%
% Returns matrix centroids where each row is an (x,y) coordinate of the
%   mouse. 1st column = X, 2nd column = Y. Length of matrix =
%   number of frames in movie
%
% Notes: Analyzes movie in chunks of chunk_size=1000 frames (smaller memory load)
%
% Updated 2015-03-10 Fori Wang

    % Default settings
    display_tracking = 0;
    trial_aware = 0;
    thresh = 20; % assumes that black mouse has RGB values <20
    
    % Check varargins
    if ~isempty(varargin)
        num_vararg = length(varargin);
        for k = 1:num_vararg
            switch varargin{k}
                case 'displayTracking'
                    display_tracking = 1;
                    figure;
                    fprintf('displayTracking ON\n');
                case 'trialAware'
                    if strfind(varargin{k+1},'.txt')
                        trial_aware = 1;
                        trial_indices = get_trial_frame_indices(varargin{k+1});
                        trial_indices = trial_indices(:,[1 4]);
                        fprintf('trialAware ON\n');
                    else
                        error('trialAware input must be .txt file.');
                    end
                case 'threshold'
                    thresh = varargin{k+1};
                    fprintf('Threshold for mouse detection provided.\n');                    
            end
        end
    end
    
    fprintf('%s: Initializing behavior video... ', datestr(now));
    behavior_vid = VideoReader(behavior_source);
    num_frames = behavior_vid.NumberofFrames;
    fprintf('Done!\n');

    % initialize centroids
    centroids = zeros(num_frames,2);
    
    % initialize variables for trial_aware mode
    lost_centroids=0;
    reassign_prev_centroid=0;
    
    % setup chunks of frames to read in movie
    chunk_size = 1000;
    frame_chunks = make_frame_chunks(num_frames,chunk_size);
    
    % setup background image: average of N frames
    num_frames_for_bg = 5000;
    fprintf('%s: Reading %d frames for background image computation... ',...
            datestr(now), num_frames_for_bg);
    bg_vid = read(behavior_vid, [1 num_frames_for_bg]);
    bg_vid = squeeze(bg_vid(:,:,1,:)); % 3D movie stack
    bg_image = mean(bg_vid,3);
    bg_image = uint8(bg_image);
    [height, width] = size(bg_image);
    fprintf('Done!\n');
    
    
    % Show the background image, and ask the user to define a polygonal
    % ROI over the main part of the maze
    h_bg = imagesc(bg_image);
    axis image;
    colormap gray;
    title(sprintf('%s: Background image',...
                  strrep(behavior_source, '_', '\_')));
    fprintf('%s: Please draw a polygonal ROI over the maze.\n', datestr(now));
    
    h_poly = impoly;
    mask_xy = getPosition(h_poly);
    pixels_to_omit = ~poly2mask(mask_xy(:,1), mask_xy(:,2), height, width);
    
    masked_bg = bg_image;
    masked_bg(pixels_to_omit) = 0;
    set(h_bg, 'CData', masked_bg);
    input('  Showing tracking ROI. Press enter to proceed >> ');
    
    % Convert the logical "off" pixels into max white in uint8, i.e. so
    % that they will not be included in the detection of a black object.
    pixels_to_omit = 255*uint8(pixels_to_omit);
    
    c_old = [0 0];
    
    for idx = 1:length(frame_chunks)

        fprintf('%s: Processing frames %d to %d...\n',...
            datestr(now), frame_chunks(idx,1), frame_chunks(idx,2));

        % read in all of the frames for the trial at the beginning
        frame_range = frame_chunks(idx,:);
        video = read(behavior_vid,frame_range);
        video = squeeze(video(:,:,1,:)); % 3D movie stack
        
        if display_tracking
            
            % plot original image on the left
            subplot(121);
            image = video(:,:,1);
            h = imagesc(image);
            title(sprintf('Original: Frames %d - %d',...
                          frame_chunks(idx,1), frame_chunks(idx,2)));
            axis image; colormap gray;
            hold on
            
            % plot initial centroids on original video
            k = plot(0,0,'k*');

            % plot initial tracking image on the right
            subplot(122);
            j = imagesc(image); colormap gray;
            title(sprintf('Tracker: Frames %d - %d',...
                          frame_chunks(idx,1),frame_chunks(idx,2)));
            hold on
            axis image;

            % plot centroids on tracking video
            l = plot(0,0,'k*');
        end

        for frame_idx = 1:size(video,3)
           
            image = video(:,:,frame_idx);
            image = max(image, pixels_to_omit); % Masked regions are replaced with white
            
            if display_tracking
                set(h,'CData',image);
                pause(0.00001);
            end
           
            % Find the mouse blob using findMouse helper function
            [final_image,s] = findMouse(bg_image,image,thresh);
            area_vector = [s.Area];
            length_vector = [s.MinorAxisLength];
            
            if display_tracking
                % Update tracking
                set(j,'CData',final_image);
                pause(0.00001);
            end

            % Save centroids data and update plot
            true_frame_idx = frame_chunks(idx,1)+frame_idx-1;           
                
            if isempty(area_vector) %sometimes blob disappears, go to previous centroid
                
                % trial_aware mode
                % if image is the first image of a trial, use next centroid
                if trial_aware && any(true_frame_idx==trial_indices(:,1))
                    reassign_prev_centroid = 1;
                    lost_centroids = lost_centroids+1;
                    centroid_color = 'm';
                elseif reassign_prev_centroid
                    % in cases mouse also lost in frames immediately
                    % following first frame
                    lost_centroids = lost_centroids+1;
                    centroid_color = 'm';
                else % use previous centroid
                    centroids(true_frame_idx,:) = c_old;
                    centroid_color = 'r';
                end
                
            else % new centroid
                [~, id] = max(length_vector); %assume mouse is the fattest blob
                c_new = s(id(1)).Centroid(1:2);
                centroids(true_frame_idx,:) = c_new;
                
                % trial_aware mode
                if reassign_prev_centroid % assign c_new to previous centroid where blob was lost
                    centroids(true_frame_idx-lost_centroids:true_frame_idx-1,:) = c_new;
                    reassign_prev_centroid = 0;
                    lost_centroids = 0;
                end
                
                c_old = c_new;
                centroid_color = 'b';
            end
            
            if display_tracking
                % Update subplots
                set(k,'XData',c_old(1),'YData',c_old(2),'color',centroid_color);
                set(l,'XData',c_old(1),'YData',c_old(2),'color',centroid_color);
                pause(0.00001);
            end
        end    
    end
    
    centroids(1,:)=centroids(2,:); %get rid of c_old [0 0] start

    % Save centroids into output text file, with "xy" extension
    [~, stem, ~] = fileparts(behavior_source);
    savename = sprintf('%s.xy', stem);
    fprintf('%s: Saving tracking results to "%s"... ', datestr(now), savename);
    fid = fopen(savename, 'w');
    for k = 1:size(centroids,1)
        fprintf(fid, '%.1f %.1f\n', centroids(k,1), centroids(k,2));
    end
    fclose(fid);
    fprintf('Done!\n');
    
end % get_mouse_XY_pos

