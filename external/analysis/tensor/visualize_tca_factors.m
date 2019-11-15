function visualize_tca_factors(tca, info, trial_meta)
% Todo:
% - Automatic identification of behavioral metadata that best separates
%   trial vector weights for a given TCA factor.

R = length(tca.lambda); % Rank of the TCA factorization
[neuron_vs, factor_assignments, cluster_boundaries] = soft_cluster_neurons(tca.U{1});
time_vs = tca.U{2};
trial_vs = tca.U{3};

num_neurons = size(neuron_vs,1);
num_trials = size(trial_vs, 1);

neuron_yrange = [0 max(neuron_vs(:))];
time_yrange = [0 max(time_vs(:))];
trial_yrange = [0 max(trial_vs(:))];

unique_days = unique(trial_meta.day);
num_unique_days = length(unique_days);

% Within-trial "time" (as opposed to index). Actual interpretation depends
% on the info.timewarp_method for trial alignment.
switch (info.timewarp_method)
    case 'align_to'
        t = info.align.axis; % Units of frames
        t0 = t(t == 0);
        t = t / 10; % Convert to seconds, assuming FPS = 10
        % The following assumes info.align.idx == 3, i.e. alignment point
        % is gate close.
        time_label = 'Time relative to gate close (s)';
        
    otherwise
        t = 1:size(time_vs,1);
        time_label = 'Time index';
end

% Prepare figure for TCA factor visualization
hf = figure; %#ok<NASGU>
% set(hf, 'MenuBar', 'none');
[ha, ~] = tight_subplot(R,3,0.01,0.01,0.01);

% For drawing vertical boundaries between different neuron clusters
num_dividers = length(cluster_boundaries);
cluster_boundaries_x = [repmat(cluster_boundaries,1,2) NaN(num_dividers,1)]';
cluster_boundaries_x = cluster_boundaries_x(:)-0.5;
cluster_boundaries_y = [repmat(neuron_yrange,num_dividers,1) NaN(num_dividers,1)]';
cluster_boundaries_y = cluster_boundaries_y(:);

all_neurons = 1:num_neurons;

for r = 1:R
    neuron_v = neuron_vs(:,r);
    time_v = time_vs(:,r);
    trial_v = trial_vs(:,r);    

    % Neuron dimension
    axes(ha((r-1)*3+1)); %#ok<*LAXES>
    neurons_in_factor = (factor_assignments==r);
    h_bar = bar(all_neurons/num_neurons, neuron_v, 'FaceColor', 'flat',...
        'EdgeColor', 'none', 'BarWidth', 1);
    h_bar.CData = repmat([0 0 0], num_neurons, 1); % Set all to black
    hold on;
    plot(cluster_boundaries_x/num_neurons, cluster_boundaries_y, 'k:');
    hold off;
    ylabel(sprintf('r = %d', r));
    
    % Time dimension
    axes(ha((r-1)*3+2));
    plot(t,time_v,'k');
    hold on;
    if (info.timewarp_method == 'align_to')
        plot(t0*[1 1], time_yrange', 'k:');
    end
    hold off;

    % Trial dimension
    axes(ha((r-1)*3+3));
    plot(trial_v, '.', 'Color', 0.4*[1 1 1], 'HitTest', 'off');
    ylabel('none');
    h_factor_trial_axes = gca;
    h_factor_trial_axes.UserData = struct(...
        'factor_idx', r,...
        'trial_coloring', 'none',...
        'x', (1:num_trials)',...
        'y', trial_v,...
        'yrange', trial_yrange,...
        'h_bar', h_bar,...
        'neurons_in_factor', neurons_in_factor);
    set(h_factor_trial_axes, 'buttonDownFcn',...
        @(h,e) recolor_trial_vector(h,e,trial_meta));
end

% Subplot formatting
%------------------------------------------------------------
all_subplots = 1:(3*R);
neuron_col_subplots = 1:3:(3*R);
time_col_subplots = 2:3:(3*R);
trial_col_subplots = 3:3:(3*R);
bottom_row_subplots = all_subplots(end-2:end);

% All plots within a column have the same scale
set(ha(neuron_col_subplots), 'XLim', [0 1], 'YLim', neuron_yrange);
set(ha(time_col_subplots), 'XLim', t([1 end]), 'YLim', time_yrange);
set(ha(trial_col_subplots), 'XLim', [1 num_trials], 'YLim', trial_yrange);

% Remove ticks and labels except for the bottom row
set(ha, 'YTickLabel', [], 'YTick', [], 'FontSize', 8);
set(ha(setdiff(all_subplots, bottom_row_subplots)), 'XTickLabel', [], 'XTick', []);

% Note: After rounding to 2 decimals, it's possible that we get ticks that
% are redundant (typically occurs with clusters with a few assigned
% neurons).
xticks = unique(round(cluster_boundaries/num_neurons, 2));
set(ha(bottom_row_subplots(1)), 'XTick', xticks,...
    'XTickLabelRotation', 90);

% Other formatting
set(ha, 'FontName', 'Arial');
set(ha(trial_col_subplots), 'YAxisLocation', 'right');

% Labels at the very bottom row
axes(ha(bottom_row_subplots(1)));
xlabel(sprintf('Neurons (%d total)', num_neurons));
axes(ha(bottom_row_subplots(2)));
xlabel(time_label);
axes(ha(bottom_row_subplots(3)));

if num_unique_days == 1
    daystr = 'on 1 day';
else
    daystr = sprintf('over %d days', num_unique_days);
end
xlabel(sprintf('Trial index (%s)', daystr));

tightfig;

end % visualize_tca_factors

function [neuron_vs, sorted_factor_assignment, cluster_boundaries] = soft_cluster_neurons(neuron_vs)
    [num_neurons, num_factors] = size(neuron_vs);
    factor_assignment = zeros(num_neurons, 1);
    
    % First, assign neurons to factors based on each neuron's weights to
    % the TCA factors.
    for k = 1:num_neurons
        neuron_factor_weights = neuron_vs(k,:);
        
        [~, sorted_inds] = sort(neuron_factor_weights, 'descend');
        factor_assignment(k) = sorted_inds(1);
    end
    [sorted_factor_assignment, sorted_neurons] = ...
        sort(factor_assignment, 'ascend');
    
    % Second, sort within each factor assignment by weight
    for l = 1:num_factors
        in_factor_l = (sorted_factor_assignment==l);
        neurons_in_factor = sorted_neurons(in_factor_l);
        neuron_factor_weights = neuron_vs(neurons_in_factor,l);
        [~, sorted_inds] = sort(neuron_factor_weights, 'descend');
        sorted_neurons(in_factor_l) = neurons_in_factor(sorted_inds);
    end
    
    neuron_vs = neuron_vs(sorted_neurons, :);
    cluster_boundaries = 1 + find(diff(sorted_factor_assignment));
end % soft_cluster_neurons

function recolor_trial_vector(h, e, trial_meta)
    data = h.UserData;
    x = data.x;
    y = data.y;
    h_bar = data.h_bar;
    neurons_in_factor = data.neurons_in_factor;
    
    if (e.Button == 3) % Right click
        new_coloring = 'none';
    else
        % Cycle through possible new coloring options
        switch (data.trial_coloring)
            case 'none'
                new_coloring = 'start';
            case 'start'
                new_coloring = 'end';
            case 'end'
                new_coloring = 'correct';
            case 'correct'
                new_coloring = 'day';
            case 'day'
                new_coloring = 'turn';
            case 'turn'
                new_coloring = 'none';
            otherwise
                new_coloring = 'none';
        end
    end
    
    % Apply new metadata coloring
    set(h, 'NextPlot', 'replacechildren');
    switch new_coloring
        case 'none'
            plot(x, y, '.', 'Color', 0.4*[1 1 1], 'HitTest', 'off');
            label_color = [0 0 0];

        case 'start'
            east_trials = strcmp(trial_meta.start, 'east');
            west_trials = strcmp(trial_meta.start, 'west');
            
            x1 = x(east_trials);
            y1 = y(east_trials);
            c1 = 'b';
            
            x2 = x(west_trials);
            y2 = y(west_trials);
            c2 = [100 150 220] / 255; % light blue
            
            plot(x1, y1, '.', 'Color', c1, 'HitTest', 'off');
            hold on;
            plot(x2, y2, '.', 'Color', c2, 'HitTest', 'off');
            hold off;
            label_color = [0 0 1];
            
        case 'end'
            north_trials = strcmp(trial_meta.end, 'north');
            south_trials = strcmp(trial_meta.end, 'south');
            dark_purple = [130 40 100] / 255;
            plot(x(north_trials), y(north_trials), '.', 'Color', dark_purple, 'HitTest', 'off');
            hold on;
            plot(x(south_trials), y(south_trials), '.', 'Color', 'm', 'HitTest', 'off');
            hold off;
            label_color = [1 0 1]; % magenta
            
        case 'correct'
            correct_trials = logical(trial_meta.correct);
            incorrect_trials = ~correct_trials;
            dark_green = [0 0.6 0];
            plot(x(correct_trials), y(correct_trials), '.', 'Color', dark_green, 'HitTest', 'off');
            hold on;
            plot(x(incorrect_trials), y(incorrect_trials), '.', 'Color', 'r', 'HitTest', 'off');
            hold off;
            label_color = dark_green;
            
        case 'day'
            unique_days = unique(trial_meta.day)'; % Row vector
            
            colors = [0.929 0.694 0.125; % Yellowish
                      0.6503 0.4858 0.0875];
            color_ind = 1;
            for unique_day = unique_days
                day_trials = find(trial_meta.day == unique_day);
                plot(x(day_trials), y(day_trials), '.', 'Color', colors(color_ind,:), 'HitTest', 'off');
                hold on;
                
                % Toggle colors
                if color_ind == 1
                    color_ind = 2;
                else
                    color_ind = 1;
                end
            end
            hold off;
            label_color = colors(1,:);
            
        case 'turn'
            left_trials = strcmp(trial_meta.turn, 'left');
            right_trials = strcmp(trial_meta.turn, 'right');
            orange = [0.85 0.3250 0.0980];
            plot(x(left_trials), y(left_trials), '.', 'Color', orange, 'HitTest', 'off');
            hold on;
            plot(x(right_trials), y(right_trials), '.', 'Color', 0.6*orange, 'HitTest', 'off');
            hold off;
            
            label_color = orange;
      
    end
    xlim(x([1 end]));
    ylim(data.yrange);
    
    if strcmp(new_coloring, 'none')
        ylabel('none', 'Color', label_color, 'FontWeight', 'normal');
    else
        ylabel(new_coloring, 'Color', label_color, 'FontWeight', 'bold');
    end
    h_bar.CData(neurons_in_factor,:) = repmat(label_color,...
        sum(neurons_in_factor), 1);
    
    % Save new state
    h.UserData.trial_coloring = new_coloring;
end