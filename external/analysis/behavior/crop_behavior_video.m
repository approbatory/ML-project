function crop_behavior_video(behavior_source)
% Spatially crop the behavior video. The cropping region is selected
%   interactively during the script. Writes the resulting file as MPEG-4
%   as '(filename)_cropped.mp4'
%
% Inputs:
%   behavior_source: Name of the behavior video.
%
% Example usage:
%   crop_behavior_video('mouse7_day07_ego-left.mp4')

behavior_video = VideoReader(behavior_source);
num_frames = behavior_video.NumberOfFrames;

% Let the user select the crop region
%------------------------------------------------------------
ref_frame = rgb2gray(read(behavior_video, 1));
imagesc(ref_frame);
colormap gray;
axis image;
title(strrep(behavior_source, '_', '\_'));
fprintf('crop_behavior_video: Please provide a rectangular region over the image.\n');
fprintf('  Double click on the rectangle when done.\n');
h_rect = imrect;
rect_params = round(wait(h_rect));

% Show the cropped image
x0 = rect_params(1); x1 = rect_params(1)+rect_params(3)-1;
y0 = rect_params(2); y1 = rect_params(2)+rect_params(4)-1;
ref_frame_cropped = ref_frame(y0:y1, x0:x1);
imagesc(ref_frame_cropped);
axis image;
title(sprintf('%s (cropped)',strrep(behavior_source, '_', '\_')));

input('crop_behavior_video: Press enter to apply crop to entire movie >> ');

% Apply crop to movie
%------------------------------------------------------------
[~, name] = fileparts(behavior_source);
output_name = sprintf('%s_cropped', name);
cropped_behavior_video = VideoWriter(output_name, 'MPEG-4');
cropped_behavior_video.Quality = 100;
cropped_behavior_video.FrameRate = 20; % FIXME: Don't hardcode

open(cropped_behavior_video);
for k = 1:num_frames
    A = read(behavior_video, k);
    A = rgb2gray(A);
    A_cr = A(y0:y1, x0:x1);
    writeVideo(cropped_behavior_video, A_cr);
    if (mod(k,1000)==0) % Show live progress
        fprintf('  %s: Frame %d of %d...\n', datestr(now), k, num_frames);
        imagesc(A_cr);
        axis image;
        title(sprintf('%s: Frame %d of %d',...
            strrep(output_name, '_', '\_'), k, num_frames));
        drawnow;
    end
end
fprintf('  %s: Done!\n', datestr(now));
close(cropped_behavior_video);