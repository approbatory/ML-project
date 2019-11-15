function M_dff = compute_dff(M, M0)
% Compute DFF of a movie with respect to a reference frame M0
%   M:  Movie to be DFF'ed [height x width x num_frames]
%   M0: Reference image [height x width]
%
% 2015 01 31 Tony Hyun Kim

num_frames = size(M,3);
if (isempty(M0))
    M0 = mean(M,3);
end

M_dff = zeros(size(M),class(M));
for k = 1:num_frames
    if (mod(k,1000)==0)
        fprintf('  Frames %d / %d completed\n', k, num_frames);
    end
    M_dff(:,:,k) = (M(:,:,k)-M0)./M0;
end