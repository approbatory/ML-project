function M_norm = norm_by_disk_filter(movie,varargin)
%Normalize every frame in the movie by a disk filtered version of itself
%
%   movie: movie matrix, [h x w x num_frames]
%   disk_radius may be passed as an input argument. Default for disk_radius
%   is 15.
% 2015 01 31 Tony Hyun Kim (Latest Revision: Hakan Inan, 15-Jan-5)
%

if isempty(varargin)
    disk_radius = 15;
elseif length(varargin) == 1
    disk_radius = varargin{1};
else
    error('Only 1 variable input argument is allowed');
end

M_norm = zeros(size(movie), 'single');
num_frames = size(movie,3);

% Apply spatial normalization
hDisk = fspecial('disk', disk_radius);

ref_idx = 1;
m_f = imfilter(single(movie(:,:,ref_idx)), hDisk, 'replicate');
m0 = mean(m_f(:));
imagesc(m_f);
xlabel('x [px]');
ylabel('y [px]');
title(sprintf('Frame %d filtered by disk of radius %d', ref_idx, disk_radius));
input('norm_by_disk_filter: Press enter to proceed >> ');

for i = 1:num_frames
    if (mod(i,1000)==0)
        fprintf('  Frames %d of %d done\n', i, num_frames);
    end
    m = single(movie(:,:,i));
    m_f = imfilter(m, hDisk, 'replicate');
    m1 = mean(m_f(:));
    m_f = m0/m1*m_f;
    M_norm(:,:,i) = m./m_f;
end
