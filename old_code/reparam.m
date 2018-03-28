function r = reparam(r)
r = 0.5-abs(0.5-r);
r = atan2(r(:,2), r(:,1))/(pi/2);
if r(end) < r(1)
    r = 1 - r;
end
end