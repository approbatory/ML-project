function draw_constant_path_analysis(ds, cell_idx)
% Used by 'view_raster_touch'

trials = ds.get_switch_trials;

pre_trials = trials.constant_pre;
post_trials = trials.constant_post;

true_trial_inds = [pre_trials; post_trials];

x_pre = 1:length(pre_trials);
x_post = x_pre(end) + (1:length(post_trials));
x_range = [1 x_post(end)];

pre = compute_features(ds, cell_idx, pre_trials);
post = compute_features(ds, cell_idx, post_trials);

% Preallocate subplots
gui.trial_times = subplot(5,2,2);
gui.mean_fluorescence = subplot(5,2,4);
gui.max_fluorescence = subplot(5,2,6);
gui.event_sum = subplot(5,2,8);
gui.event_counts = subplot(5,2,10);

% Behavioral trial times
draw_stem(gui.trial_times, pre.times, post.times);
ylabel('Trial times (s)');
title(sprintf('Constant path: Pre (blue, %d) vs. Post (red, %d)',...
              length(pre_trials), length(post_trials)));

% Fluorescence
draw_stem(gui.mean_fluorescence, pre.mean_fluorescence, post.mean_fluorescence);
ylabel('Mean fluorescence');

draw_stem(gui.max_fluorescence, pre.max_fluorescence, post.max_fluorescence);
ylabel('Max fluorescence');

% Event amplitude sum
draw_stem(gui.event_sum, pre.event_sum, post.event_sum);
ylabel('\Sigma Event amplitudes');

% Event counts
draw_stem(gui.event_counts, pre.num_events, post.num_events);
xlabel('Trial index');
ylabel('Event counts');

    function draw_stem(h_sp, pre_vals, post_vals)
        subplot(h_sp);
        stem(x_pre, pre_vals, 'b.', 'ButtonDownFcn', {@select_trial, h_sp});
        hold on;
        stem(x_post, post_vals, 'r.', 'ButtonDownFcn', {@select_trial, h_sp});
        hold off;
        xlim(x_range);
        tick_inds = x_range(1):5:x_range(end);
        xticks(tick_inds);
        xticklabels(num2cell(true_trial_inds(tick_inds)));
        ylim(compute_ylim(pre_vals, post_vals));
        grid on;
    end % draw_stem

    function select_trial(h_stem, e, h_sp)
        x = e.IntersectionPoint(1);
        y = h_stem.YData(h_stem.XData==x);
        t = true_trial_inds(x);
        
        subplot(h_sp);
        hold on;
        plot(x,y,'ko');
        text(x+0.5,y,sprintf('Trial %d', t));
        hold off;
    end

end % draw_constant_path_analysis

function y_range = compute_ylim(pre_vals, post_vals)
    all_vals = [pre_vals(:); post_vals(:)];
    m = min(min(all_vals), 0);
    M = max(all_vals);
    
    if (m ~= M)
        y_range = [m M] + 0.1*(M-m)*[-1 1];
    else
        y_range = m + [-0.5 0.5];
    end
end