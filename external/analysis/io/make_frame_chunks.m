function [ frame_chunks, num_chunks ] = make_frame_chunks( num_frames , chunk_size)
% Creates frame chunks so that get_mouse_XY_pos can go
%   through the behavior video in chunks (not enough memory to
%   load the whole video)
%
%   splits up frames into chunks of size chunk_size
%     For example, [frame_indices] = make_frame_indices(50000, 1000)
%       returns:           [1    1000
%                          1001  2000
%                          2001  3000 ...] etc. until the last frame
%
% 2015-02-26 Fori Wang

    x = (1:chunk_size:num_frames)';
    y = [x(2:end)-1; num_frames];
    frame_chunks = [x y];
    num_chunks = size(frame_chunks,1);
    
end

