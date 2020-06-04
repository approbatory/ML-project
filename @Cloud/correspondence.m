function correspondence(o)

load asymp_snr_values.mat

[mat, sem, varnames] = o.predictor_matrix;

[asymp_ratio, asymp_ratio_conf] = ...
    uncertain_divide(asymp_snr, asymp_snr_conf,...
    asymp_snr_shuf, asymp_snr_shuf_conf);

[asymp_worst, asymp_worst_conf] = ...
    worst_divide(asymp_snr, asymp_snr_conf,...
    asymp_snr_shuf, asymp_snr_shuf_conf);


conf = sem * 1.96;

for j = 1:numel(varnames)
    [pearson(j,1), pearson_p(j,1)] = corr(mat(:,j), asymp_ratio(:), 'type', 'Pearson');
    [spearman(j,1), spearman_p(j,1)] = corr(mat(:,j), asymp_ratio(:), 'type', 'Spearman');
    [kendall(j,1), kendall_p(j,1)] = corr(mat(:,j), asymp_ratio(:), 'type', 'Kendall');
    [~, adjr2(j,1)] = Utils.regress_line(mat(:,j), asymp_ratio(:));    
end

keyboard


end

function [quotient, quotient_uncertainty] = uncertain_divide(x, xc, y, yc)
quotient = x./y;
quotient_uncertainty = abs(x./y).*sqrt((xc./x).^2 + (yc./y).^2);
end

function [quotient, quotient_uncertainty] = worst_divide(x, xc, y, yc)
lx = x - xc;
ux = x + xc;
ly = y - yc;
uy = y + yc;

lq = lx./uy;
uq = ux./ly;

quotient = (uq + lq) ./ 2;
quotient_uncertainty = (uq - lq) ./ 2;
end