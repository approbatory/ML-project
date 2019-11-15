function F = compute_fluorescence_stats(M, varargin)
% Compute fluorescence stats (min, mean, max) of a movie M on a 
% frame-by-frame basis.
%
% Inputs:
%   M: [height x width x num_frames] movie to be examined
%

trim = 0;
if ~isempty(varargin)
    for k = 1:length(varargin)
        switch lower(varargin{k})
            case 'trim' % Trim borders
                trim = varargin{k+1};
        end
    end
end

if (ischar(M)) % 'M' is the name of a HDF5 file
    [movie_size, ~] = get_movie_info(M);
    num_frames = movie_size(3);
    
    F = zeros(size(M,3), 3, 'single');
    
    chunk_size = 2500;
    [frame_chunks, num_chunks] = make_frame_chunks(num_frames, chunk_size);
    
    write_idx = 0;
    for chunk_idx = 1:num_chunks
        fprintf('%s: Examining chunk %d of %d from %s...\n',...
            datestr(now), chunk_idx, num_chunks, M);
        
        M_chunk = load_movie_from_hdf5(M, frame_chunks(chunk_idx,:));
        
        for k = 1:size(M_chunk,3)
            write_idx = write_idx + 1;
            m = M_chunk((1+trim):(end-trim),...
                  (1+trim):(end-trim),...
                  k);
            m = single(m(:));
            F(write_idx,:) = [min(m) mean(m) max(m)];
        end
    end
    fprintf('%s: Done!\n', datestr(now));
else % 'M' is the movie matrix itself
    num_frames = size(M,3);
    
    F = zeros(size(M,3), 3, 'single');
    
    for frame_idx = 1:num_frames
        if (mod(frame_idx,2500)==0)
            fprintf('%s: Frames %d of %d examined...\n',...
                datestr(now), frame_idx, num_frames);
        end
        m = M((1+trim):(end-trim),...
              (1+trim):(end-trim),...
               frame_idx);
        m = single(m(:));
        F(frame_idx,:) = [min(m) mean(m) max(m)];
    end
end




