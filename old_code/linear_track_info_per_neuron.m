%% script parameters
visualize = true; %use to make intermediate plots

%% loading file

data_options = { 'Mouse-2024-20150315_090450-linear-track-TracesAndEvents.mat',...
    'Mouse-2024-20150317_202708-linear-track-TracesAndEvents.mat'};
chosen_option = data_options{2};

path_to_data = fullfile('../linear_track', chosen_option);

data_struct = load(path_to_data);
tracesEvents = data_struct.tracesEvents;

%% constructing dataset
full_X = tracesEvents.rawProb;%tracesEvents.rawTraces;
full_y = tracesEvents.position;

[forward_mask, backward_mask, cm_per_pix] = select_directions(full_y);
%take forward trials
my_X = full_X(forward_mask,:);
my_y = full_y(forward_mask,:);

%% divide into bins

fullLen = 118;
midLen = range(my_y(:,1)) .* cm_per_pix; %cm
num_bins = 20;
my_binner = @(y) gen_place_bins(y, num_bins, midLen);
[my_ks, my_centers] = my_binner(my_y);

%%
if visualize
    figure;
    subplot(2,2,1);
    imagesc(full_X.');
    xlabel('Frames, \Deltat=50ms')
    ylabel('Neurons');
    title('Full neural traces');
    subplot(2,2,2);
    imagesc(my_X.');
    xlabel('Frames, \Deltat=50ms')
    ylabel('Neurons');
    title('Forward trials'' neural traces');
    
    subplot(2,2,3);
    y_min = min(full_y(:,1)); range_y = range(full_y(:,1));
    renorm_y = @(y) (y-y_min)./range_y.*fullLen;
    coord_to_plot = renorm_y(full_y(:,1));
    plot(coord_to_plot);
    hold on;
    partial_y = coord_to_plot;
    partial_y(~forward_mask) = nan;
    plot(partial_y, 'r');
    legend full 'forward trials'
    xlabel('Frames, \Deltat=50ms');
    ylabel('Track coordinate (cm)');
    title('Full track coordinate');
    subplot(2,2,4);
    coord_to_plot = renorm_y(my_y(:,1));
    plot(coord_to_plot, '-k');
    hold on;
    plot(coord_to_plot, '.b');
    partial_bins_y = coord_to_plot;
    partial_bins_y(mod(my_ks,2) == 0) = nan;
    plot(partial_bins_y, '.r');
    for b_n = 1:20
        refline(0,b_n*midLen./num_bins);
    end
    for f_n = 1:10
        line([f_n f_n].*length(coord_to_plot)./10, ylim);
    end
    legend 'position trace' 'even bins' 'odd bins'
    xlabel('Frames, \Deltat=50ms');
    ylabel('Track coordinate (cm)');
    title('Forward trials'' track coordinate');
    
    figure;
    [sorted_my_ks, bin_order] = sort(my_ks);
    normified_my_X = my_X - min(my_X);
    normified_my_X = normified_my_X ./ sum(normified_my_X);
    neurons_mean_ks = normified_my_X.' * my_ks;
    [~, neurons_ks_order] = sort(neurons_mean_ks);
    sorted_my_X = my_X(bin_order,neurons_ks_order);
    %sorted_my_X = (sorted_my_X - mean(sorted_my_X))./std(sorted_my_X);
    sorted_my_X = (sorted_my_X - min(sorted_my_X))./sum((sorted_my_X - min(sorted_my_X)));
    imagesc(sorted_my_X.');
    hold on;
    divisions = find(diff(sorted_my_ks)~=0);
    for d = divisions.'
        line([d d], ylim, 'Color', [1 1 1]);
    end
    title('Neural responses per bin (min set to 0 and normalized)');
    xlabel('Frames, divided into bins');
    ylabel('Neurons, ascending bin tuning');
end


%% selections of cells

delta_neurons = 10;
total_neurons = size(my_X,2);
neuron_nums = delta_neurons:delta_neurons:total_neurons;
if neuron_nums(end) ~= total_neurons
    neuron_nums = [neuron_nums total_neurons];
end
if delta_neurons > 1
    neuron_nums = [1 neuron_nums];
end

num_samples = 1;%20; %how many times to run each number of neurons
tot_ticker = tic;
%for c_ix = numel(neuron_nums):-1:1
res_size = [numel(neuron_nums) num_samples];
empty_field = cell(res_size);
res = struct('MSE',empty_field,'mean_err',empty_field,'MSE_train',empty_field,'mean_err_train',empty_field,'num_cells',empty_field, 'ks', empty_field, 'pred', empty_field);
res_shuf = res;
%parfor c_ix = 1:numel(neuron_nums)
for c_ix = numel(neuron_nums):-1:1
    number_of_neurons = neuron_nums(c_ix);
    sub_ticker = tic;
    for rep_ix = 1:num_samples
        cell_subset_X = my_X(:, randperm(total_neurons) <= number_of_neurons);
        cell_subset_X_shuf = shuffle(cell_subset_X, my_ks);
        res(c_ix, rep_ix) = decoding_reporter(cell_subset_X, my_ks, my_centers);
        res_shuf(c_ix, rep_ix) = decoding_reporter(cell_subset_X_shuf, my_ks, my_centers);
        fprintf('%d cells :- %d ', number_of_neurons, rep_ix);
    end
    toc(sub_ticker);
    fprintf('\n %d cells done (ix=%d), err=%.2f +- %.2f cm || err_shuf=%.2f +- %.2f cm\n', number_of_neurons, c_ix, mean(mean(cat(1,res(c_ix, :).mean_err),2)), std(mean(cat(1,res(c_ix, :).mean_err),2))./sqrt(num_samples), mean(mean(cat(1,res_shuf(c_ix, :).mean_err),2)), std(mean(cat(1,res_shuf(c_ix, :).mean_err),2))./sqrt(num_samples));
end
toc(tot_ticker);
%% plotting
if visualize
    for c_ix = numel(neuron_nums):-1:1
        err_per_div =  mean(cat(1,res(c_ix, :).mean_err),2);
        err_mean(c_ix) = mean(err_per_div);
        err_errb(c_ix) = std(err_per_div)./sqrt(length(err_per_div));
        
        err_per_div_shuf = mean(cat(1,res_shuf(c_ix, :).mean_err),2);
        err_mean_shuf(c_ix) = mean(err_per_div_shuf);
        err_errb_shuf(c_ix) = std(err_per_div_shuf)./sqrt(length(err_per_div_shuf));
        
    end
    figure;
    errorbar(neuron_nums, err_mean, err_errb);
    hold on;
    errorbar(neuron_nums, err_mean_shuf, err_errb_shuf);
    legend unshuffled shuffled
    xlabel('Number of cells');
    ylabel('Mean error (cm)');
    title('Mean decoding error vs. Number of cells used');
end

%% saving
if ~exist('records', 'dir')
    mkdir('records');
end
save(sprintf('records/lin_track_2024_0317_%s.mat', timestring), 'res', 'res_shuf');
%% decoding function
function res = decoding_reporter(X, ks, centers)
k_fold = 10;
alg = my_algs('ecoclin');
%alg = my_algs('lda');

train_slices = ceil((1:length(ks))./length(ks).*k_fold) ~= (1:k_fold).';

MSE = @(k,p) mean(sum((centers(k,:)-centers(p,:)).^2,2));
mean_err = @(k,p) mean(sum(abs(centers(k,:)-centers(p,:)),2));

res.num_cells = size(X,2);
res.ks = zeros(size(ks));
res.pred = zeros(size(ks));

for i_fold = 1:size(train_slices,1)
    train_slice = train_slices(i_fold,:);
    test_slice = ~train_slice;
    X_train = X(train_slice,:); X_test = X(test_slice,:);
    ks_train = ks(train_slice); ks_test = ks(test_slice);
    model = alg.train(X_train, ks_train);
    pred_train = alg.test(model, X_train);
    pred_test = alg.test(model, X_test);
    res.MSE(i_fold) = MSE(ks_test, pred_test);
    res.mean_err(i_fold) = mean_err(ks_test, pred_test);
    res.ks(test_slice) = ks_test;
    res.pred(test_slice) = pred_test;
    
    res.MSE_train(i_fold) = MSE(ks_train, pred_train);
    res.mean_err_train(i_fold) = mean_err(ks_train, pred_train);
end

end