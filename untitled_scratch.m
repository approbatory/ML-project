DecodeTensor.aggregate_results('save_dir', 'records_confidence',...
    'table_name', 'pairwise', 'db_file', 'pairwise.db',...
    'field_names', {'Mouse', 'Setting', 'BinA', 'BinB', 'CorrectConfidence'},...
    'create_command',...
    'create table pairwise(Mouse text, Setting text, BinA int, BinB int, CorrectConfidence real);',...
    'save_varname', 'records');


%%
command_template = @(isbwA, isbwB, offset, setting, mouse)...
            sprintf(['select avg(CorrectConfidence) from pairwise'...
            ' where BinA between %d*20+1 and %d*20+20 and'...
            ' BinB between %d*20+1 and %d*20+20 and'...
            ' BinB = BinA + %d*20 + %d and'...
            ' Setting = ''%s'' and Mouse = ''%s'' group by BinA, BinB;'],...
            isbwA, isbwA, isbwB, isbwB, isbwB && ~isbwA, offset, setting, mouse);
        
      
%%

[source_path, mouse_name] = DecodeTensor.default_datasets(3);
l_ = load(source_path);
[fast_frames, fast_regular_frames, bins, dirbins] = ...
    DecodeTensor.aux_sel(l_.tracesEvents.position, DecodeTensor.default_opt);

num_reg_frames = sum(fast_regular_frames);
X = l_.tracesEvents.rawTraces(fast_regular_frames, :);
ks = dirbins(fast_regular_frames);

X_irreg = l_.tracesEvents.rawTraces(fast_frames & ~fast_regular_frames, :);
ks_irreg = dirbins(fast_frames & ~fast_regular_frames);

tr_subset = randperm(num_reg_frames) <= floor(num_reg_frames/2);
X_tr = X( tr_subset, :);
X_te = X(~tr_subset, :);
ks_tr = ks( tr_subset);
ks_te = ks(~tr_subset);

alg = my_algs('ecoclin');
model = alg.train(X_tr, ks_tr);
preds_reg = alg.test(model, X_te);
preds_irreg = alg.test(model, X_irreg);

binsize = 118/20;
mean_err_func = @(ks, ps) mean(abs(ceil(ks/2) - ceil(ps/2))) * binsize;

reg_mean_err = mean_err_func(ks_te, preds_reg);
irreg_mean_err = mean_err_func(ks_irreg, preds_irreg);

fprintf('Train on 1/2 regular, test on other 1/2 regular:\t%f cm error\nTrain on 1/2 regular, test on irregular:\t%f cm error\n', reg_mean_err, irreg_mean_err);

%% pooled plot, normed at max neurons
dbfile = 'decoding.db';
conn = sqlite(dbfile); %remember to close it
mouse_list = {'Mouse2010','Mouse2012',...
    'Mouse2023','Mouse2026',...
    'Mouse2019','Mouse2028',...
    'Mouse2024','Mouse2022'};
%% mid
for i = 1:numel(mouse_list)
    [nP{i}, mP{i}, eP{i}] = ...
        DecodingPlotGenerator.get_errors('NumNeurons', conn,...
        mouse_list{i}, 'unshuffled', 'IMSE', 'max');
    [nPs{i}, mPs{i}, ePs{i}] = ...
        DecodingPlotGenerator.get_errors('NumNeurons', conn,...
        mouse_list{i}, 'shuffled', 'IMSE', 'max');
    
    index_of_180 = find(nPs{i} == 180,1);
    norm_by = mPs{i}(index_of_180);
    mPs_normed{i} = mPs{i} ./ norm_by;
    ePs_normed{i} = ePs{i} ./ norm_by;
    
    %index_of_180 = find(nPs{i} == 180,1);
    %norm_by = mPs{i}(index_of_180);
    mP_normed{i} = mP{i} ./ norm_by;
    eP_normed{i} = eP{i} ./ norm_by;
end
%%
lookup_at_value = @(val, val_col, res_col) cell2mat(...
    cellfun(@(n,m) m(n == val), val_col, res_col,...
    'UniformOutput', false).');
%%
m_at_n = @(n) mean(lookup_at_value(n, nP, mP_normed));
e_at_n = @(n) sqrt(var(lookup_at_value(n, nP, mP_normed)) +...
    mean(lookup_at_value(n, nP, eP_normed).^2));

m_at_n_s = @(n) mean(lookup_at_value(n, nPs, mPs_normed));
e_at_n_s = @(n) sqrt(var(lookup_at_value(n, nPs, mPs_normed)) +...
    mean(lookup_at_value(n, nPs, ePs_normed).^2));

%%
figure; hold on;
nn = [1 (30:30:500)];
mm = arrayfun(m_at_n, nn);
mm_s = arrayfun(m_at_n_s, nn);
ee = arrayfun(e_at_n, nn);
ee_s = arrayfun(e_at_n_s, nn);
errorbar(nn, mm, ee, 'b');
errorbar(nn, mm_s, ee_s, 'r');
%%
figure; hold on;
for i = 1:numel(mouse_list)
    plot(nP{i}, mP_normed{i}, 'b');
end
for i = 1:numel(mouse_list)
    plot(nPs{i}, mPs_normed{i}, 'r');
end
%% closing section
conn.close;