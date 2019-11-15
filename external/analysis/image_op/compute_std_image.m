function [S, A] = compute_std_image(M)
% Compute the per-pixel standard deviation in a way that is more memory
% efficient than S = std(M,[],3);

[height, width, num_frames] = size(M);
A = zeros(height, width);
A2 = zeros(height, width);

for k = 1:num_frames
    A = A + M(:,:,k);
    A2 = A2 + M(:,:,k).^2;
end

A = A / num_frames;
A2 = A2 / num_frames;
S = sqrt(A2-A.^2);