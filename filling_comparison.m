%comparing filling perf
PL = 16;


ds = quick_ds('../open_field/om3-0327', 'nocells');

%fillings = {'binary', 'box', 'copy_zeroed', 'copy', 'traces'};
fillings = {'box', 'traces'};

%algs = my_algs({'ecoclin_onevsall'}, {'original', 'shuf'}, true, 0);
algs{1} = [my_algs('mvnb2'), my_algs('ecoclin_onevsall'), my_algs('ecoclin_onevsall', 'shuf')];
algs{2} = [my_algs('lda'), my_algs('gnb'), my_algs('ecoclin_onevsall'), my_algs('ecoclin_onevsall', 'shuf')];


for ix = 1:numel(fillings)
    [X, ks, errf] = ds_dataset(ds, 'filling', fillings{ix},...
        'openfield', true, 'sparsify', ~strcmp(fillings{ix}, 'traces'));
    for i = 1:numel(algs{ix})
        [train_err{ix,i}, test_err{ix,i}] = evaluate_alg(algs{ix}(i),...
            X, ks, 'eval_f', errf, 'train_frac', 0.7, 'par_loops', PL);
    end
    fprintf('Completed filling: %s\n', fillings{ix});
end

%%
fillings{3} = 'copy (zeroed)';
bin_cm = 5.12;
conv_to_cm = @(Q) cellfun(@(x) x.*5.12 ,Q, 'UniformOutput', false);

cm_train_err = conv_to_cm(train_err);
cm_test_err = conv_to_cm(test_err);
figure;
berr(fillings, cm_test_err, cm_train_err, {algs.name});

%plot labels:
title('Preprocessing Effect on Place Decoding');
xlabel('Trace preprocessing method');
ylabel('Mean bin distance, scaled by bin size (cm)');


dirname = 'graphs/place_decoding/open_field_hpc/filling_comparison';
fname = ['filling_comparison_' timestring '.png'];

set(gcf,'PaperUnits','inches','PaperPosition',2*[0 0 4 3]);
print(fullfile(dirname,fname), '-dpng', '-r100');





%% FILLING_COMPARISON V2 : using new binning
PL = 4;

ds = quick_ds('../open_field/om3-0327', 'nocells');
fillings = {'box', 'traces'};
algs{1} = [my_algs('mvnb2'), my_algs('ecoclin_onevsall'), my_algs('ecoclin_onevsall', 'shuf')];
algs{2} = [my_algs('lda'), my_algs('gnb'), my_algs('ecoclin_onevsall'), my_algs('ecoclin_onevsall', 'shuf')];
for ix = 1:numel(fillings)
    [X, ~, ~] = ds_dataset(ds, 'filling', fillings{ix},...
        'openfield', true, 'sparsify', ~strcmp(fillings{ix}, 'traces'));
    y = ds.trials.centroids;
    %[ks, centers, dims, scXY, ~] = gen_place_bins(ds, 8, 46);
    %errf_xy = @(xy,p) sqrt(mean(sum((xy - centers(p(:),:)).^2,2)));
    for i = 1:numel(algs{ix})
    %    [train_err{ix,i}, test_err{ix,i}] = evaluate_alg(algs{ix}(i),...
    %        X, ks, 'eval_f_xy', errf_xy, 'XY', scXY, 'train_frac', 0.7, 'par_loops', PL);
        [meas_train{ix,i}, meas_test{ix,i}, ~] = evala(algs{ix}(i), X, y, @(y) gen_place_bins(y,8,46),...
            'repeats', PL, 'verbose', true, 'split', 'nonlocal');
    end
    fprintf('Completed filling: %s\n', fillings{ix});
end


figure;
ylab = 'RMS error (cm)';
bsep({'Place Decoders'' Performance on Box (nonlocal split)', 'Place Decoders'' Performance on Traces (nonlocal split)'}, meas_test, meas_train, Inf, {{algs{1}.short},{algs{2}.short}}, ylab);
%% FILLING_COMPARISON V3 : using bin center interpolation -- DEBUGGING
PL = 16;
PAR = true;
%ds = quick_ds('../open_field/om2-0405', 'nocells');
%ds = quick_ds('../open_field/c11m1d25', 'nocells');
fillings = {'traces'};
algs = [my_algs('lda'), my_algs('dlda'), my_algs('dqda')];
for ix = 1:numel(fillings)
    %[X, ~, ~] = ds_dataset(ds, 'filling', fillings{ix},...
    %    'openfield', true, 'sparsify', ~strcmp(fillings{ix}, 'traces'));
    %y = ds.trials.centroids;
    y = E_T{6}.tracesEvents.position; %(T, 2)
    if max(max(y)) > 1e3
        continue;
    end
    X = E_T{6}.tracesEvents.rawTraces; %(T, C)
    %[ks, centers, dims, scXY, ~] = gen_place_bins(ds, 8, 46);
    %errf_xy = @(xy,p) sqrt(mean(sum((xy - centers(p(:),:)).^2,2)));
    for i = 1:numel(algs)
    %    [train_err{ix,i}, test_err{ix,i}] = evaluate_alg(algs{ix}(i),...
    %        X, ks, 'eval_f_xy', errf_xy, 'XY', scXY, 'train_frac', 0.7, 'par_loops', PL);
        [meas_train{ix,i}, meas_test{ix,i}, ~] = evala(algs(i), X, y, @(y) gen_place_bins(y,20,120),...
            'repeats', PL, 'verbose', true, 'interp', false, 'split', 'nonlocal', 'use_par', PAR);
        [meas_train_interp{ix,i}, meas_test_interp{ix,i}, ~] = evala(algs(i), X, y, @(y) gen_place_bins(y,20,120),...
            'repeats', PL, 'verbose', true, 'interp', true, 'split', 'nonlocal', 'use_par', PAR);
    end
    fprintf('Completed filling: %s\n', fillings{ix});
end

% 
% figure;
% ylab = 'RMS error (cm)';
% bsep({'Place Decoders'' Performance on Traces'}, meas_test, meas_train, Inf, {{algs.short}}, ylab);
% 
% 
% figure;
% ylab = 'RMS error (cm)';
% bsep({'Place Decoders'' Performance on Traces, with bin interpolation'}, meas_test_interp, meas_train_interp, Inf, {{algs.short}}, ylab);


figure;
ylab = 'RMS error (cm)';
bsep({'Place Decoders'' Performance on Traces', 'Place Decoders'' Performance on Traces, with bin interpolation'},...
    [meas_test;meas_test_interp], [meas_train;meas_train_interp], 8, {{algs.short},{algs.short}}, ylab);
