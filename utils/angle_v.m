function [theta, cos_angle] = angle_v(a,b)
a = a(:); b = b(:);
cos_angle = (a.' * b) ./ norm(a) ./ norm(b);
theta = acosd(cos_angle);
end