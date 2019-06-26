d = DecodeTensor(6);

[X, ks] = d.get_dataset;
progressbar('temp_corrs');
for i = 1:20
    cv_model = fitcecoc(X, ks, 'CrossVal', 'on', 'KFold', 2);
    ks_hat = kfoldPredict(cv_model);
    
    
    mean_err(i) = mean(abs(ceil(ks/2) - ceil(ks_hat/2))) .* 5.9;
    
    fprintf('%d:: Mean error is: %.3f\n', i, mean_err(i));
    progressbar(i/20);
end

fprintf('Total:: Mean error: %.3f +- %.3f\n', mean(mean_err), std(mean_err)./sqrt(20));

%%
o = Analyzer('../linear_track/Mouse2022/Mouse-2022-20150326-linear-track/Mouse-2022-20150326_093722-linear-track-TracesAndEvents.mat');
opt = DecodeTensor.default_opt;
[~,~,tr_start,tr_end,~,~,ks] = DecodeTensor.new_sel(o.data.y.raw.full, opt);

tr_mask = zeros(numel(ks),1);
tr_lens = zeros(1,numel(tr_start));
for tr_i = 1:numel(tr_start)
    tr_mask(tr_start(tr_i):tr_end(tr_i)) = tr_i;
    tr_lens(tr_i) = tr_end(tr_i) - tr_start(tr_i) + 1;
end

ks_cut = ks(tr_mask~=0);
X_cut = o.data.X.full(tr_mask~=0, :);

%%
progressbar('bad temcorrs');
for i = 1:20
    cv_model_bad = fitcecoc(X_cut, ks_cut, 'CrossVal', 'on', 'KFold', 2);
    ks_cut_hat = kfoldPredict(cv_model_bad);
    
    mean_err_bad(i) = mean(abs(ceil(ks_cut/2) - ceil(ks_cut_hat/2))) .* 5.9;
    fprintf('%d:: Mean error is: %.3f\n', i, mean_err_bad(i));
    progressbar(i/20);
end

fprintf('Total:: Mean error: %.3f +- %.3f\n', mean(mean_err_bad), std(mean_err_bad)./sqrt(20));
