%% define save file
save_file = '../linear_track/pablo_data_ds.mat';
output_file = '../linear_track/pablo_data_output.mat';
%load the data into a ds array
loaded_res = load(save_file);
my_ds = loaded_res.my_ds;

%% relevant variables about each day
number_of_cells = [my_ds.num_cells];
number_of_frames = [my_ds.full_num_frames];

number_of_traverses = zeros(1,numel(my_ds));
for ix = 1:numel(my_ds)
    [~, bw] = select_directions(my_ds(ix));
    number_of_traverses(ix) = sum(diff(bw) == -1);
end


%% decoding outcomes containers
%multiple dimensions
% dim 1 : samples
% dim 2 : test/train
% dim 3 : decoder
% dim 4 : day index
% dim 5 : forwards/backwards
% dim 6 : normal/shuffle
algs = [my_algs({'gnb', 'lda'}),my_algs('ecoclin', {'shuf', 'original'})'];
num_samples = 16;

res.dims = {'samples', 'test/train', 'decoder', 'day_index', 'forwards/backwards', 'normal/shuffle'};

res.dims_meaning = ...
    {1:num_samples, {'test', 'train'}, {algs.short}, 1:numel(my_ds),...
    {'forwards', 'backwards'}, {'normal', 'shuffle'}};
res.size = cellfun(@length, res.dims_meaning);
res.array = zeros(res.size);
%% COPYING OLD RESULTS
for ix = 1:res.size(4)
    for i = 1:res.size(3)
        if length(tr_err{ix,i}) == num_samples
            res.array(:,2,i,ix,1,1) = tr_err{ix,i};
        end
        if length(te_err{ix,i}) == num_samples
            res.array(:,1,i,ix,1,1) = te_err{ix,i};
        end
        
        if length(tr_err_sh{ix,i}) == num_samples
            res.array(:,2,i,ix,1,2) = tr_err_sh{ix,i};
        end
        if length(te_err_sh{ix,i}) == num_samples
            res.array(:,1,i,ix,1,2) = te_err_sh{ix,i};
        end
        
        
        if length(tr_err_back{ix,i}) == num_samples
            res.array(:,2,i,ix,2,1) = tr_err_back{ix,i};
        end
        if length(te_err_back{ix,i}) == num_samples
            res.array(:,1,i,ix,2,1) = te_err_back{ix,i};
        end
        
        if length(tr_err_sh_back{ix,i}) == num_samples
            res.array(:,2,i,ix,2,2) = tr_err_sh_back{ix,i};
        end
        if length(te_err_sh_back{ix,i}) == num_samples
            res.array(:,1,i,ix,2,2) = te_err_sh_back{ix,i};
        end
    end
end
%% new results
midLen = 0.8*120; %cm
num_bins = 20;
for norm_shuf_6 = 1:res.size(6)
    for for_back_5 = 1:res.size(5)
        for day_idx_4 = 1:res.size(4)
            for alg_idx_3 = 1:res.size(3)
                my_alg = algs(alg_idx_3);
                my_X = my_ds(day_idx_4).trials.traces.';
                my_y = my_ds(day_idx_4).trials.centroids;
                [sel_fw, sel_bw] = select_directions(my_y);
                if for_back_5 == 1
                    my_X = my_X(sel_fw,:);
                    my_y = my_y(sel_fw,:);
                elseif for_back_5 == 2
                    my_X = my_X(sel_bw,:);
                    my_y = my_y(sel_bw,:);
                else
                    error('wrong option');
                end
                my_binner = @(y) gen_place_bins(y, num_bins, midLen);
                
                my_ticker = tic;
                [res.array(:,2,alg_idx_3,day_idx_4,for_back_5,norm_shuf_6),...
                res.array(:,1,alg_idx_3,day_idx_4,for_back_5,norm_shuf_6)] = ...
                    evala(my_alg, my_X, my_y, my_binner, 'split', 'nonlocal',...
                    'repeats', res.size(1), 'verbose', true,...
                    'shufboth', norm_shuf_6==2, 'use_par', true, 'errfunc', 'mean_dist');
                time_taken = toc(my_ticker);
                fprintf(': : %d %d %d %d\t%.2f s\n', alg_idx_3,day_idx_4,for_back_5,norm_shuf_6, time_taken);
            end
        end
    end
end
%% saving output file
save(output_file, '-struct', 'res');
%% plotting
mean_errs = mean(res.array);
errb_errs = std(res.array)./sqrt(res.size(1));

figure;
subplot(2,1,1);
errnbar(squeeze(mean_errs(1,1,:,:,1,1)).', squeeze(errb_errs(1,1,:,:,1)).');
hold on;
errnbar(squeeze(mean_errs(1,2,:,:,1,1)).', squeeze(errb_errs(1,2,:,:,1)).', true);
legend(algs.short);
xlabel('Linear track sessions (mouse 2028)');
ylabel('mean error (cm)');
title('Place decoding from traces (forward)');
set(gca, 'XTickLabel', res.dims_meaning{4});
set(gca, 'XTick', 1:res.size(4));
ylim([0 35]);

subplot(2,1,2);
errnbar(squeeze(mean_errs(1,1,:,:,2,1)).', squeeze(errb_errs(1,1,:,:,1)).');
hold on;
errnbar(squeeze(mean_errs(1,2,:,:,2,1)).', squeeze(errb_errs(1,2,:,:,1)).', true);
legend(algs.short);
xlabel('Linear track sessions (mouse 2028)');
ylabel('mean error (cm)');
title('Place decoding from traces (backward)');
set(gca, 'XTickLabel', res.dims_meaning{4});
set(gca, 'XTick', 1:res.size(4));
ylim([0 35]);

suptitle('Linear track place decoding, separately for forward/backward passes');
%% DIAGNOSTIC: DIRECTION SELECTION
seldec = @select_directions;
ix = 18;
[fw, bw] = seldec(my_ds(ix));
figure;
x_coord = my_ds(ix).trials.centroids(:,1);
plot(x_coord, 'b');
hold on
fw = double(fw); bw = double(bw);
fw(fw == 0) = nan;
bw(bw == 0) = nan;
plot(x_coord.*fw, 'r');
plot(x_coord.*bw, 'g');
title('Selecting directions');
xlabel('Frames');
ylabel('X coordinate');
legend coordinate forward backward

%%
function [fw, bw] = select_directions(XY)
if isstruct(XY)
    XY = XY.trials.centroids;
end
vthresh = 15 / 20; %~3cm/s
track_coord = XY(:,1);
velocity = diff(track_coord);
smooth_velocity = medfilt1(velocity, 21);
is_forward = smooth_velocity > vthresh;
is_backward = smooth_velocity < -vthresh;

track_range = range(track_coord);
track_min = min(track_coord);
track_max = max(track_coord);
bottom_tenth = track_min + track_range/10;
top_tenth = track_max - track_range/10;

in_between = (track_coord > bottom_tenth) & (track_coord < top_tenth);
fw = [false; is_forward] & in_between;
bw = [false; is_backward] & in_between;
end

%% OLD FUNC
% function [fw, bw] = select_directions(XY)
% if isstruct(XY)
%     XY = XY.trials.centroids;
% end
% track_coord = XY(:,1);
% velocity = diff(track_coord);
% smooth_velocity = smooth(velocity, 200);
% is_forward = smooth_velocity > 0;
% 
% track_range = range(track_coord);
% track_min = min(track_coord);
% track_max = max(track_coord);
% bottom_tenth = track_min + track_range/10;
% top_tenth = track_max - track_range/10;
% 
% in_between = (track_coord > bottom_tenth) & (track_coord < top_tenth);
% fw = [false; is_forward] & in_between;
% bw = [false; (~is_forward)] & in_between;
% end