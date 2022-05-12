function f5e
%% Figure 5e
% Depends on files: default_store.mat
% Raw numbers for:
% - Bar2021X, Bar2021Y, Bar2026X, Bar2026Y
rng(0);

data = place_field_widths;
make_xlsx(data, 'f5e');
end

%% helper functions
function data = place_field_widths
org = Org;

filt_21 = strcmp(org.mouse, 'Mouse2021');
filt_26 = strcmp(org.mouse, 'Mouse2026');

mus_all = org.fetch('mus');
mus_21 = mus_all(filt_21);
mus_26 = mus_all(filt_26);

hwhm_21 = cell2mat(cellfun(@calc_hwhm, mus_21, 'UniformOutput', false));
hwhm_26 = cell2mat(cellfun(@calc_hwhm, mus_26, 'UniformOutput', false));

fwhm_21 = 2*5.9*hwhm_21;
fwhm_26 = 2*5.9*hwhm_26;

figure;
histogram(fwhm_21, 'BinWidth', 2*5.9, 'Normalization', 'probability');

[data.Bar2021Y, data.Bar2021X] = ...
    histcounts(fwhm_21, 'BinWidth', 2*5.9, 'Normalization', 'probability');

xlabel 'PF widths [cm]'
ylabel 'prob. density'
title 'Mouse2021'
xlim([0 120]);
ylim([0 0.33]);

figure;
histogram(fwhm_26, 'BinWidth', 2*5.9, 'Normalization', 'probability');

[data.Bar2026Y, data.Bar2026X] = ...
    histcounts(fwhm_26, 'BinWidth', 2*5.9, 'Normalization', 'probability');

xlabel 'PF widths [cm]'
ylabel 'prob. density'
title 'Mouse2026'
xlim([0 120]);
ylim([0 0.33]);
end

function hwhm = calc_hwhm(mus)
mus = cell2mat(mus);
[n_cells, n_bins] = size(mus);
assert(n_bins == 40);
[hwhm_r, hwhm_l] = deal(zeros(1, n_cells));
for i = 1:n_cells
    mu_r = mus(i, 1:20);
    mu_l = mus(i, 21:40);
    
    if all(mu_r == mu_r(1))
        hwhm_r(i) = nan;
    else
        hwhm_r(i) = one_hwhm(mu_r);
    end
    
    if all(mu_l == mu_l(1))
        hwhm_l(i) = nan;
    else
        hwhm_l(i) = one_hwhm(mu_l);
    end
end
hwhm = [hwhm_r hwhm_l];
hwhm = hwhm(~isnan(hwhm));
end

function hwhm = one_hwhm(x)
x = x - min(x);
[mx, ix_max] = max(x);
hmx = mx/2;

pre_x = x(1:ix_max);
post_x = x(ix_max:end);

left_bound = find(pre_x < hmx, 1, "last");
if isempty(left_bound)
    left_hw = -inf;
else
    left_hw = ix_max - left_bound;
end

right_bound = find(post_x < hmx, 1, "first");
if isempty(right_bound)
    right_hw = -inf;
else
    right_hw = right_bound - 1;
end

hwhm = max(left_hw, right_hw);
if hwhm < 0
    assert(all(x >= hmx));
    warning('Value never dipped below half max.');
    hwhm = round(numel(x)/2);
end
end