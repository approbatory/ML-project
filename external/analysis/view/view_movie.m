function view_movie(M, varargin)
% Displays the frames of a movie matrix M [height x row x num_frames]
%   (Note: also works with a single image)

movie_clim = [];
use_mask = 0;
use_outline = 0;
poi = [];
num_repeats = inf;
rescale_each_frame = false;
for k = 1:length(varargin)
    vararg = varargin{k};
    if ischar(vararg)
        switch lower(vararg)
            case {'repeat', 'repeats'}
                num_repeats = varargin{k+1};
            case 'rescale'
                rescale_each_frame = true;
            case 'mask'
                % Pixels with logical 1 in the provided mask will be displayed
                use_mask = true;
                ds = varargin{k+1};
                mask = ds.get_mask;
            case 'poi'
                poi = varargin{k+1};
            case 'clim'
                movie_clim = varargin{k+1};
            case 'boundary'
                use_outline = 1;
                ds = varargin{k+1};
        end
    end
end

num_frames = size(M,3);

if ~isempty(movie_clim) % If CLim is provided, use it
    h = imagesc(M(:,:,1), movie_clim);
else
    if (rescale_each_frame || isa(M, 'uint16'))
        % Raw movies (e.g. uint16) are rescaled by default
        h = imagesc(M(:,:,1));
        mask_val = 0; % FIXME
    else % Otherwise, use common CLim scaling
        movie_clim = compute_movie_scale(M);
        mask_val = movie_clim(1);
        h = imagesc(M(:,:,1), movie_clim);
    end
end

axis image;
truesize;
colormap gray;
xlabel('x [px]');
ylabel('y [px]');

if use_outline
    hold on;
    cell_indices = find(ds.is_cell);
    for cell_idx = cell_indices
        boundary = ds.cells(cell_idx).boundary;
        plot(boundary(:,1), boundary(:,2), 'g');
    end
    hold off;
end

if ~isempty(poi)
    hold on;
    plot(poi(:,1), poi(:,2), 'm*');
    hold off;
end

num_playbacks = 1;
while (num_playbacks <= num_repeats) 
    for k = 1:num_frames
        title(sprintf('Frame %d of %d', k, num_frames));
        m = M(:,:,k);
        if use_mask
            m(mask) = mask_val;
        end
        set(h, 'CData', m);
        drawnow;
    end
    num_playbacks = num_playbacks + 1;
end
