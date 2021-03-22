function [robust_std_jitter, deviations] = estimate_position_jitter(id, plotting)
if nargin == 1
    plotting = false;
end

m = SessManager;
d = m.cons_usable(id, true);
TE = load(d{1});
TE = TE.tracesEvents;

[cpp, vel, trial_start, trial_end, trial_direction, track_bins, track_dir_bins] = ...
    DecodeTensor.new_sel(TE.position, DecodeTensor.default_opt);
pos = TE.position(:,1);

deviations = [];

bin_s = 7; bin_e = 13;
bin_s_a = 2*bin_s-1; bin_s_b = 2*bin_e;
bin_e_a = 2*bin_e-1; bin_e_b = 2*bin_s;

for i = 1:numel(trial_start)

    tr_s = trial_start(i);
    tr_e = trial_end(i);
    
    tr_bins = track_dir_bins(tr_s:tr_e);
    tr_pos = cpp * pos(tr_s:tr_e);
    
    lin_s = find((tr_bins == bin_s_a) | (tr_bins == bin_s_b), 1);
    lin_e = find((tr_bins == bin_e_a) | (tr_bins == bin_e_b), 1, 'last');
    
    if isempty(lin_s) || isempty(lin_e)
        continue;
    end
    
    
    frame_idx = lin_s:lin_e;
    real_pos = tr_pos(frame_idx);
    
    training_subset = mod(frame_idx, 2) == 0;
    testing_subset = ~training_subset;
    
    if sum(training_subset) < 4
        fprintf('Skipping trial with only %d training samples', sum(training_subset));
        continue;
    end
    
    linfit = Utils.regress_line(frame_idx(training_subset),...
        real_pos(training_subset));
    fit_pos = linfit(frame_idx);
    %fit_pos = linspace(real_pos(1), real_pos(end), numel(real_pos))';
    
    if plotting
        figure;
        subplot(1,2,1);
        plot(real_pos, 'o'); hold on; plot(fit_pos, 'o'); legend real fit
        title(sprintf('Trial %d', i));
        subplot(1,2,2);
        plot(real_pos - fit_pos, 'o');
        pause
    end
    
    devs = real_pos - fit_pos;
    %devs = devs(2:end-1);
    devs = devs(testing_subset);
    deviations = [deviations; devs];
end

%rms_jitter = sqrt(mean(deviations.^2));
robust_std_jitter = iqr(deviations) / (2*sqrt(2)*erfinv(1/2));