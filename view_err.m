function view_err(ds, poss, err, err_map, label, varargin)
%VIEW_ERR Summary of this function goes here
%   Detailed explanation goes here
save_figs = ~isempty(varargin) && strcmp(varargin{1}, 'save');
POS_LABEL = 'position along arm (0.25=middle of arm)';
figure;
plot(poss, err, '-x');
title(sprintf('%s: End arm decoding (NB) vs. position', label));
xlabel(POS_LABEL);
ylabel('Leave-1-out x-val error');
if save_figs
    print(sprintf('Multinomial_NB_error_curve_for_%s.png',label), '-dpng');
end

figure;
imagesc(1 - err_map - (1-err_map).*0.1.*mod((1:size(err_map,1))',2));
colormap('gray');
[n_trials, n_pos] = size(err_map);
for i = 1:n_trials
    if ds.trials(i).correct
        c = 'g';
    else
        c = 'r';
    end
    rectangle('Position', [(n_pos-0.5) (i-0.5) 2 1], 'FaceColor', c);
end

ticks = get(gca, 'XTick') + 1;
set(gca, 'XTick', ticks, 'XTickLabel', poss(ticks));
t = text(n_pos+2, n_trials/4, 'non-rewarded in red');
t.Rotation = -90;
xlabel(POS_LABEL);
ylabel('Trial number');
title(sprintf('%s: Map of errors (black)', label));
if save_figs
    print(sprintf('Multinomial_NB_error_map_for_%s.png',label), '-dpng');
end

end

