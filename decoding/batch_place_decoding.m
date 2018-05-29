function res = batch_place_decoding(dayset_labels, varargin)
p = inputParser;

p.addRequired('dayset_labels', @(x) iscell(x) || ischar(x));
p.addParameter('filling', 'traces', @ischar);
p.addParameter('include_spoof', true, @islogical);
%p.addParameter('include_unmoving', true, @islogical);
p.addParameter('parloops', 1, @isnumeric);
p.addParameter('train_frac', 0.7, @isnumeric);
p.addParameter('msg', 'No message', @ischar);
p.addParameter('pickup', [], @(x) ischar(x) || isempty(x));
p.parse(dayset_labels, varargin{:});
dayset_labels = p.Results.dayset_labels;
filling = p.Results.filling;
include_spoof = p.Results.include_spoof;
%include_unmoving = p.Results.include_unmoving;
parloops = p.Results.parloops;
train_frac = p.Results.train_frac;
msg = p.Results.msg;
pickup = p.Results.pickup;

if strcmp(filling, 'traces')
    algs = my_algs({'gnb', 'ecoclin_onevsall'}, {'original', 'shuf'}, true, 0);
else
    algs = my_algs({'mvnb2', 'ecoclin_onevsall'}, {'original', 'shuf'}, true, 0);
end

if isempty(pickup)
    pickup = timestring;
end
save_filename = sprintf('batch_place_decoding_res_%s.mat', pickup);
if ~exist(save_filename, 'file')
    is_new = true;
    save(save_filename, 'msg', 'dayset_labels');
else
    S = load(save_filename, 'res');
    res = S.res;
    is_new = false;
end

daysets = auto_dayset(dayset_labels);
if is_new
    for j = 1:numel(daysets)
        res{j}.finished = false(numel(daysets{j}), numel(algs));
    end         
end
for j = 1:numel(daysets)
    %dayset_label = dayset_labels{j};
    process_dayset(daysets{j}, j);
end


    function process_dayset(dayset, j)
        res{j}.dayset = dayset;
        res{j}.algs = algs;
        for ix = 1:numel(dayset)
            disp(dayset(ix).day);
            d = dayset(ix);
            ds = quick_ds(fullfile(d.directory, d.day), {'deprobe', 'nocells'});
            % ds.num_trials == 1
            [X, ks, errf] = ds_dataset(ds, 'filling', filling,...
                'openfield', ds.num_trials == 1,...
                'sparsify', ~strcmp(filling, 'traces'));
            spoof_samp_gen = @() shuffle(X, ks);
            subset = subset_moving(ds);
            for i = 1:numel(algs)
                if res{j}.finished(ix,i)
                    fprintf('SKIPPING: at indices j%d ix%d i%d\n', j, ix, i);
                    continue;
                end
                [res{j}.train_err{ix,i}, res{j}.test_err{ix,i},...
                    res{j}.sub_train_err{ix,i}, res{j}.sub_test_err{ix,i},...
                    res{j}.moving_train_err{ix,i}, res{j}.moving_test_err{ix,i}] = ...
                    evaluate_alg(algs(i), X, ks, 'eval_f', errf,...
                    'subset', subset, 'train_frac', train_frac,...
                    'par_loops', parloops);
                if include_spoof
                    [res{j}.spoof_train_err{ix,i}, res{j}.spoof_test_err{ix,i},...
                        res{j}.spoof_sub_train_err{ix,i}, res{j}.spoof_sub_test_err{ix,i},...
                        res{j}.spoof_moving_train_err{ix,i}, res{j}.spoof_moving_test_err{ix,i}] = ...
                        evaluate_alg(algs(i), spoof_samp_gen, ks, 'eval_f', errf,...
                        'subset', subset, 'train_frac', train_frac,...
                        'par_loops', parloops, 'X_is_func', true);
                end
                res{j}.finished(ix,i) = true;
                save(save_filename, '-append', 'res', 'msg', 'dayset_labels');
                fprintf('Saved: at indices j%d ix%d i%d\n', j, ix, i);
            end
        end
    end
end