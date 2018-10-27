function [theta, cos_angle] = angle_v(a,b, use_abs)
if ~exist('use_abs', 'var')
    use_abs = false;
end
a = a(:); b = b(:);
cos_angle = (a.' * b) ./ norm(a) ./ norm(b);
if use_abs
    cos_angle = abs(cos_angle);
end
theta = acosd(cos_angle);
end