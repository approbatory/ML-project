function view_err(ds, poss, err, err_map, label, varargin)
%VIEW_ERR Plot error curves / error map from the output of decode_end_nb
%optionally save the plots, or hide the appearance of the plots in windows
save_figs = false;
show_figs = true;
savedir = '.';
for k = 1:length(varargin)
    if ischar(varargin{k})
        switch varargin{k}
            case 'save'
                save_figs = true;
                savedir = varargin{k+1};
            case 'hide'
                show_figs = false;
        end
    end
end
POS_LABEL = 'position along arm (0.25=middle of arm)';
if show_figs
    figure;
else
    figure('Visible','Off');
end
plot(poss, err, '-x');
title(sprintf('%s: End arm decoding (NB) vs. position', label));
xlabel(POS_LABEL);
ylabel('Leave-1-out x-val error');
if save_figs
    fname = fullfile(savedir, sprintf('Multinomial_NB_error_curve_for_%s.png',label));
    print(fname, '-dpng');
end

if show_figs
    figure;
else
    figure('Visible','Off');
end
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
    fname = fullfile(savedir,sprintf('Multinomial_NB_error_map_for_%s.png',label));
    print(fname, '-dpng');
end

end

