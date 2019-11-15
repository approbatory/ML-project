function traces = extract_traces(movie_source,rec_source)
% Extract traces from movie source and rec source. The inputs can be either
% the paths to the sources or the sources themselves. Specifically, movie
% can be fed as the path, 3D movie matrix or the 2D (pixels flattened)
% movie matrix. Likewise, the rec source might be the path to the rec file
% or the 'filters' structure in the rec file (or its transpose).

regul = 0; % regularization term for LS reconstruction

% Determine if movie source is the link to source or the movie matrix
if ischar(movie_source)
    M = load_movie(movie_source);
    [h,w,t] = size(M);
    num_pixels = h*w;
    M = reshape(M,h*w,t);
else
    if ismatrix(movie_source)
        M = movie_source;
        num_pixels = size(M,1);
    else
         [h,w,t] = size(M);
         M = reshape(M,h*w,t);
         num_pixels = h*w;
    end
end

% Determine if rec source is the link to source or the actual filters
if ischar(rec_source)
    rec = load(rec_source);
    filters = rec.filters;
    num_filters = size(filters,3);
    F = reshape(filters,h*w,num_filters);
else
    if ismatrix(rec_source)
        [d1,d2] = size(rec_source);
        if d1==num_pixels
            F = rec_source;
        elseif d2==num_pixels
            F = rec_source';
        else
            error('Sizes of movie matrix and filters must be consistent.')
        end
    else
        num_filters = size(rec_source,3);
        F = reshape(rec_source,h*w,num_filters);
    end
end

idx_nonzero = find(sum(F,2)>0);
M_small = M(idx_nonzero,:);
F_small = F(idx_nonzero,:);
T = (F_small'*F_small+regul*eye(size(F,2))) \  (F_small'*M_small);
traces = T';
