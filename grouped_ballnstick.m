function [p,h] = grouped_ballnstick(labels, ys, errs, varargin)
%n = size(labels, 2);

averaged_ys = cellfun(@weighted_avg, ys, errs);
errs_w = cellfun(@weighted_errs, ys, errs);
[p,h] = ballnstick(labels{1}, labels{2},...
    averaged_ys(1,:), averaged_ys(2,:),...
    errs_w(1,:), errs_w(2,:), varargin{:});


end

function res = weighted_avg(ys, errs)
w = 1./errs;
w = w./sum(w);
res = sum(ys.*w);
end

function res = weighted_errs(ys, errs)
w = 1./errs;
w = w./sum(w);
mu = weighted_avg(ys, errs);
v1 = sum(w);
v2 = sum(w.^2);
res = sqrt(sum(w.*(ys - mu).^2)./(v1 - v2./v1));
end