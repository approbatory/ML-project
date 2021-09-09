function [quotient, quotient_uncertainty] = uncertain_divide(x, xc, y, yc)
quotient = x./y;
quotient_uncertainty = abs(x./y).*sqrt((xc./x).^2 + (yc./y).^2);
end