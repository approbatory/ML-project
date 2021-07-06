function res = bin_integration_time_checker(i_usable, n_divs_array, integration_frames_array, reps, alg)
if ~exist('alg', 'var')
    alg = my_algs('ecoclin');
end

ntype = 'events_transients';

t_ = tic;

[coordwise_error, binwise_error] = deal(zeros(numel(n_divs_array), numel(integration_frames_array), reps));


%i_usable = 69;
if isstruct(i_usable)
    tracesEvents = i_usable;
else
    sm = SessManager;
    data_file = sm.data_path(i_usable, 'Usable');
    load(data_file, 'tracesEvents');
end
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

NUM_NDIVS = numel(n_divs_array);
NUM_INT_FRAMES = numel(integration_frames_array);
%n_divs = 5;
%progressbar('#Reps', '#Bins', '#Integration frames');
parfor r_ix = 1:reps
    for n_ix = 1:NUM_NDIVS
        n_divs = n_divs_array(n_ix);
        milestones = linspace(0,opt.total_length, n_divs+1);
        %milestones = opt.total_length./n_divs .* (0.5:1:n_divs);
        key_events = zeros(numel(trial_start), n_divs); %frame numbers
        for i = 1:numel(trial_start)
            st = trial_start(i);
            en = trial_end(i);
            pos_segment = cm_position(st:en);
            for j = 1:n_divs
                %ms = milestones(j);
                %dist_to_ms = abs(pos_segment - ms);
                %[~, closest_frame] = min(dist_to_ms);
                %key_events(i,j) = st + closest_frame - 1;
                ms_begin = milestones(j); ms_end = milestones(j+1);
                frames_within = pos_segment >= ms_begin & pos_segment < ms_end;
                if ~any(frames_within)
                    key_events(i,j) = nan;
                else
                    num_frames_within = sum(frames_within);
                    random_frame = find(frames_within,1) + randi(num_frames_within) - 1;
                    key_events(i,j) = st + random_frame - 1;
                end
            end
        end
        
        %integration_frames = 4;
        for i_ix = 1:NUM_INT_FRAMES
            integration_frames = integration_frames_array(i_ix);
            N = size(tracesEvents.(ntype),2);
            K = n_divs;
            T = numel(trial_start);
            neural_features = zeros(N,K,T);
            [position_target, bin_target] = deal(zeros(K,T));
            selector = (1:integration_frames) - floor(integration_frames/2) - 1;
            max_frame = size(tracesEvents.(ntype),1);
            for k_i = 1:K
                for t_i = 1:T
                    frame = key_events(t_i, k_i);
                    if isnan(frame)
                        neural_features(:,k_i, t_i) = nan;
                        position_target(k_i, t_i) = nan;
                        bin_target(k_i, t_i) = nan;
                    else
                        frame_range = frame + selector;
                        frame_range = min(max_frame, frame_range);
                        frame_range = max(1, frame_range);
                        %ind = sub2ind([K T], k_i, t_i);
                        neural_features(:,k_i, t_i) = mean(tracesEvents.(ntype)(frame_range,:),1);
                        position_target(k_i, t_i) = cm_position(frame);
                        bin_target(k_i, t_i) = 2*k_i - (tr_dir(t_i)==1);
                    end
                end
            end
            
            ideal_bin_position = zeros(2*K,1);
            for k_val = 1:2*K
                x_ = position_target(bin_target==k_val);
                ideal_bin_position(k_val) = nanmean(x_(:));
            end
            
            %%
            %[~,~,ps,ks,~,~]
            %mean_err = DecodeTensor.decode_tensor(neural_features, tr_dir, opt.total_length/n_divs, my_algs('ecoclin'), false, [], []);
            %fprintf('The mean error was %f cm\n', mean_err);
            %TODO use the analog position for error calculation
            %%
            [T1, d1, T2, d2, division] = DecodeTensor.holdout_half(neural_features, tr_dir);
            tr_ix_1 = find( division);
            tr_ix_2 = find(~division);
            pos_1 = position_target(:, division);
            pos_2 = position_target(:,~division);
            
            [sup_X1, sup_ks1, f_x1] = DecodeTensor.tensor2dataset(T1, d1);
            sup_pos1 = pos_1(:);
            %[sup_X1, f_x1] = denan_rows(sup_X1);
            %sup_ks1 = sup_ks1(f_x1);
            sup_pos1 = sup_pos1(f_x1);
            
            [sup_X2, sup_ks2, f_x2] = DecodeTensor.tensor2dataset(T2, d2);
            sup_pos2 = pos_2(:);
            %[sup_X2, f_x2] = denan_rows(sup_X2);
            %sup_ks2 = sup_ks2(f_x2);
            sup_pos2 = sup_pos2(f_x2);
            
            binsize = opt.total_length/n_divs;
            binwise_err_func = @(ks, ps) sqrt(mean((ceil(ks/2) - ceil(ps/2)).^2) * binsize.^2);%mean(abs(ceil(ks/2) - ceil(ps/2))) * binsize;
            coordwise_err_func = @(pos, ideal_pred_pos) sqrt(mean((pos - ideal_pred_pos).^2));%mean(abs(pos - ideal_pred_pos));
            
            %alg = my_algs('ecoclin');
            
            %train on part 1
            model = alg.train(sup_X1, sup_ks1);
            sup_ps2 = alg.test(model, sup_X2);
            binwise_error2 = binwise_err_func(sup_ks2, sup_ps2);
            coordwise_error2 = coordwise_err_func(sup_pos2, ideal_bin_position(sup_ps2));
            
            if any(isnan(coordwise_error2(:)))
                keyboard;
            end
            
            
            %train on part 2
            model = alg.train(sup_X2, sup_ks2);
            sup_ps1 = alg.test(model, sup_X1);
            binwise_error1 = binwise_err_func(sup_ks1, sup_ps1);
            coordwise_error1 = coordwise_err_func(sup_pos1, ideal_bin_position(sup_ps1));
            
            if any(isnan(coordwise_error1(:)))
                keyboard;
            end
            
            binwise_error(n_ix, i_ix, r_ix) = mean([binwise_error1 binwise_error2]);
            coordwise_error(n_ix, i_ix, r_ix) = mean([coordwise_error1 coordwise_error2]);
            
            fprintf('ndivs=%d frames=%d reps=%d\nCoordwise RMS error:\t%.2f cm\nBinwise RMS error:\t%.2f cm\n\n',...
                n_divs, integration_frames, r_ix, coordwise_error(n_ix, i_ix, r_ix), binwise_error(n_ix, i_ix, r_ix));
            %progressbar([], [], i_ix / numel(integration_frames_array));
        end
        %progressbar([], n_ix / numel(n_divs_array), []);
    end
    %progressbar(r_ix/reps, [], []);
end

res.c_err = coordwise_error;
res.b_err = binwise_error;
res.nbins = n_divs_array;
res.nframes = integration_frames_array;
res.reps = reps;
res.usable_index = i_usable;


time_taken = toc(t_);
fprintf('Took %d s to complete\n', time_taken);

x_ = mean(coordwise_error,3);
if isvector(x_) || isscalar(x_)
    return;
end


figure;
%[NN, II] = meshgrid(n_divs_array, integration_frames_array);
surf(integration_frames_array, n_divs_array, mean(coordwise_error,3));
xlabel 'Integration frames'
ylabel 'Number of bins'
zlabel 'Mean error (cm)'
title 'Coordwise error'


figure;
%[NN, II] = meshgrid(n_divs_array, integration_frames_array);
surf(integration_frames_array, n_divs_array, mean(binwise_error,3));
xlabel 'Integration frames'
ylabel 'Number of bins'
zlabel 'Mean error (cm)'
title 'Binwise error'



end

function [x, f] = denan_rows(x)
    f = ~any(isnan(x),2);
    x = x(f,:);
end