sm = SessManager;
data_file = sm.data_path(69, 'Usable');
load(data_file, 'tracesEvents');
opt = DecodeTensor.default_opt;
[cpp, ~, trial_start, trial_end, tr_dir, track_bins, track_dir_bins] = ...
    DecodeTensor.new_sel(tracesEvents.position, opt);
%%

track_coord = tracesEvents.position(:,1);
track_begin = prctile(track_coord, opt.cutoff_p);
track_end = prctile(track_coord, 100 - opt.cutoff_p);

cm_position = (track_coord - track_begin) .* cpp;
cm_position(cm_position < 0) = 0;
s_time = (1:numel(cm_position)) / 20;

n_divs = 5;
%milestones = linspace(0,opt.total_length, n_divs);
milestones = opt.total_length./n_divs .* (0.5:1:n_divs);
key_events = zeros(numel(trial_start), numel(milestones)); %frame numbers
for i = 1:numel(trial_start)
    st = trial_start(i);
    en = trial_end(i);
    pos_segment = cm_position(st:en);
    for j = 1:numel(milestones)
        ms = milestones(j);
        dist_to_ms = abs(pos_segment - ms);
        [~, closest_frame] = min(dist_to_ms);
        key_events(i,j) = st + closest_frame - 1;
    end
end

integration_frames = 4;

N = size(tracesEvents.rawTraces,2);
K = n_divs;
T = numel(trial_start);
neural_features = zeros(N,K,T);
[position_target, bin_target] = deal(zeros(1,K,T));
selector = (1:integration_frames) - floor(integration_frames/2) - 1;
max_frame = size(tracesEvents.rawTraces,1);
for k_i = 1:K
    for t_i = 1:T
        frame = key_events(t_i, k_i);
        frame_range = frame + selector;
        frame_range = min(max_frame, frame_range);
        frame_range = max(1, frame_range);
        %ind = sub2ind([K T], k_i, t_i);
        neural_features(:,k_i, t_i) = mean(tracesEvents.rawTraces(frame_range,:),1);
        position_target(k_i, t_i) = cm_position(frame);
        bin_target(k_i, t_i) = k_i + K*(tr_dir(t_i)==-1);
    end
end

ideal_bin_position = zeros(1,2*K);
for k_val = 1:2*K
    x_ = position_target(bin_target==k_val);
    ideal_bin_position(k_val) = mean(x_(:));
end

%%
mean_err = DecodeTensor.decode_tensor(neural_features, tr_dir, opt.total_length/n_divs, my_algs('ecoclin'), false, N, floor((T-3)/2));
fprintf('The mean error was %f cm\n', mean_err);
%TODO use the analog position for error calculation