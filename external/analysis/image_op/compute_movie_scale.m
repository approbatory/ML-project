function clim = compute_movie_scale(M)
% Compute an appropriate viewing scale (CLim) for the provided movie

[height, width, ~] = size(M);
maxVec = reshape(max(M,[],3), height*width, 1);
minVec = reshape(min(M,[],3), height*width, 1);
quantsMax = quantile(maxVec,[0.85,0.87,0.9,0.93,0.95]);
quantsMin = quantile(minVec,[0.85,0.87,0.9,0.93,0.95]);

clim = [mean(quantsMin) mean(quantsMax)];
clim_range = clim(2)-clim(1);
clim = clim + 0.1*clim_range*[-1 1];
