function A = compute_mean_image(M, frame_inds)
% Compute the mean image over specified 'frame_inds' of the movie 'M'. Is
% far more memory efficient that A = mean(M(:,:,frame_inds),3);

[height, width, num_frames] = size(M);
A = zeros(height, width);

if ~exist('frame_inds', 'var')
    frame_inds = 1:num_frames;
end

if iscolumn(frame_inds)
    frame_inds = frame_inds'; % Want row
end

for k = frame_inds
    A = A + M(:,:,k);
end

A = A / length(frame_inds);