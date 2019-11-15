function compare_registered_movie(M_orig, M_reg, varargin)
% Simultaneously display original (unregistered) and registered movies,
%   with linked XY ranges. Optional input will repeat the movie.
%
% Movies are [h x w x num_frames]. The original and registered movies
%   are expected to have the same dimensions.
%
% 2015 01 30 Tony Hyun Kim
if isempty(varargin)
    num_repeats = 20;
else
    num_repeats = varargin{1};
end

num_frames = size(M_orig, 3);
M_scale = [min(M_orig(:)) max(M_orig(:))];

ax_orig = subplot(1,2,1);
h_orig = imagesc(M_orig(:,:,1), M_scale);
axis image;
colormap gray;
title('Original');

ax_reg = subplot(1,2,2);
h_reg = imagesc(M_reg(:,:,1), M_scale);
axis image;
colormap gray;
title('Registered');

linkaxes([ax_orig, ax_reg], 'xy');

for r = 1:num_repeats
    for k = 1:num_frames
        set(h_orig, 'CData', M_orig(:,:,k));
        set(h_reg, 'CData', M_reg(:,:,k));
        drawnow;
    end
end
