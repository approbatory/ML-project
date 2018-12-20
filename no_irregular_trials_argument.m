reg_mean_err = zeros(1, 8); irreg_mean_err = reg_mean_err; combined_mean_err = reg_mean_err;
mouse_names = cell(1,8);
for dispatch_index = 1:8
    [source_path, mouse_name] = DecodeTensor.default_datasets(dispatch_index);
    mouse_names{dispatch_index} = mouse_name;
    l_ = load(source_path);
    if strcmp(mouse_name, 'Mouse2022')
        track_coord = l_.tracesEvents.position(91:end, :);
        X_full = l_.tracesEvents.rawTraces(91:end, :);
    else
        track_coord = l_.tracesEvents.position;
        X_full = l_.tracesEvents.rawTraces;
    end
    [fast_frames, fast_regular_frames, bins, dirbins] = ...
        DecodeTensor.aux_sel(track_coord, DecodeTensor.default_opt);
    
    %num_reg_frames = sum(fast_regular_frames);
    X = X_full(fast_regular_frames, :);
    [num_reg_frames, total_neurons] = size(X);
    ks = dirbins(fast_regular_frames);
    
    X_irreg = X_full(fast_frames & ~fast_regular_frames, :);
    ks_irreg = dirbins(fast_frames & ~fast_regular_frames);
    
    
    while false
        neuron_series = [1 (50:50:total_neurons) total_neurons];
        reg_mean_err = zeros(20, numel(neuron_series));
        irreg_mean_err = zeros(20, numel(neuron_series));
        for n = 1:numel(neuron_series)
            parfor r = 1:20 %reps
                [reg_mean_err(r,n), irreg_mean_err(r,n)] =...
                    irreg_decode(X, ks, X_irreg, ks_irreg, neuron_series(n));
                fprintf(['n_neu=%d, rep=%d\nTrain on 1/2 regular,'...
                    ' test on other 1/2 regular:\t%f cm error\n'...
                    'Train on 1/2 regular,'...
                    ' test on irregular:\t\t%f cm error\n\n'],...
                    neuron_series(n), r, ...
                    reg_mean_err(r,n), irreg_mean_err(r,n));
            end
        end
        
        save no_irregular_trials_argument.mat reg_mean_err irreg_mean_err source_path mouse_name
    end
    %% doing it for all mice
    d = dispatch_index;
    r = 1;
    [reg_mean_err(r, d), irreg_mean_err(r, d), combined_mean_err(r, d)] = ...
        irreg_decode(X, ks, X_irreg, ks_irreg, total_neurons);
    fprintf(['d=%d, rep=%d\nTrain on 1/2 regular,'...
        ' test on other 1/2 regular:\t%f cm error\n'...
        'Train on 1/2 regular,'...
        ' test on irregular:\t\t%f cm error\n'...
        'Train on 1/2 regular,'...
        ' test on both:\t\t\t%f cm error\n\n'],...
        d, r, ...
        reg_mean_err(r,d), irreg_mean_err(r,d), combined_mean_err(r,d));
end
%%
while false
errb = @(x)std(x)./sqrt(length(x));
reg_m = mean(reg_mean_err); reg_e = errb(reg_mean_err);
irg_m = mean(irreg_mean_err); irg_e = errb(irreg_mean_err);

figure;
handles = [];
h1 = DecodingPlotGenerator.errors_plotter(neuron_series, reg_m, reg_e, 'unshuffled', 'DisplayName', 'Regular trials'' frames tested');
h2 = DecodingPlotGenerator.errors_plotter(neuron_series, irg_m, irg_e, 'shuffled', 'DisplayName', 'Irregular trials'' frames tested');
set(gca, 'YScale', 'log');
legend([h2 h1]);
%legend Location east;
legend boxoff;
xlabel 'Number of neurons';
ylabel 'Mean decoding error (cm)';
title(sprintf('Regular vs. irregular trials, %s', mouse_name));
end

figure;
bar([reg_mean_err; irreg_mean_err; combined_mean_err]');
legend regular irregular combined

%xlabel Mouse
ylabel 'Decoding error (cm)'
title 'Regular vs. irregular frames, training on 50% of both, testing on 50% regular, 50% irregular, and combined'
set(gca, 'XTickLabel', mouse_names);
ylim([0 40]);

function [reg_mean_err, irreg_mean_err, combined_mean_err] = irreg_decode(X, ks, X_irreg, ks_irreg, n_neurons)
train_on_both = true;
[num_reg_frames, total_neurons] = size(X);
tr_subset = (1:num_reg_frames) <= floor(num_reg_frames/2);
neuron_subset = randperm(total_neurons, n_neurons);
X_tr = X( tr_subset, neuron_subset);
X_te = X(~tr_subset, neuron_subset);
ks_tr = ks( tr_subset);
ks_te = ks(~tr_subset);

num_irreg_frames = size(X_irreg,1);
tr_subset_irreg = (1:num_irreg_frames) <= floor(num_irreg_frames/2);
X_irreg_tr = X_irreg( tr_subset_irreg, neuron_subset);
X_irreg_te = X_irreg(~tr_subset_irreg, neuron_subset);
ks_irreg_tr = ks_irreg(tr_subset_irreg);
ks_irreg_te = ks_irreg(~tr_subset_irreg);

X_tr_combined = [X_tr ; X_irreg_tr];
ks_tr_combined = [ks_tr ; ks_irreg_tr];
X_te_combined = [X_te ; X_irreg_te];
ks_te_combined = [ks_te ; ks_irreg_te];

alg = my_algs('ecoclin');
if train_on_both
    model = alg.train(X_tr_combined, ks_tr_combined);
else
    model = alg.train(X_tr, ks_tr);
end
preds_reg = alg.test(model, X_te);
preds_irreg = alg.test(model, X_irreg_te);
combined_preds = alg.test(model, X_te_combined);
combined_ks = ks_te_combined;

binsize = 118/20;
mean_err_func = @(ks, ps) mean(abs(ceil(ks/2) - ceil(ps/2))) * binsize;

reg_mean_err = mean_err_func(ks_te, preds_reg);
irreg_mean_err = mean_err_func(ks_irreg_te, preds_irreg);
combined_mean_err = mean_err_func(combined_ks, combined_preds);
% fprintf(['Train on 1/2 regular,'...
%     ' test on other 1/2 regular:\t%f cm error\n'...
%     'Train on 1/2 regular,'...
%     ' test on irregular:\t%f cm error\n'],...
%     reg_mean_err, irreg_mean_err);
end