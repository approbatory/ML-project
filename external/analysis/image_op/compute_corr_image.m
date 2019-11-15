function C = compute_corr_image(M, point)

[height, width, num_frames] = size(M);
point = round(point);

x = max(1,point(1)); x = min(width,x);
y = max(1,point(2)); y = min(height,y);

ref_trace = squeeze(M(y,x,:));
ref_trace = 1/norm(ref_trace)^2 * ref_trace;

Mr = reshape(M, height*width, num_frames);
C = reshape(Mr * ref_trace, height, width);