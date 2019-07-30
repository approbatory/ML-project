function res = confusion_mat_runner(index)
tic
D_T = DecodeTensor.cons_filt(index);
n_reps = 20;
num_trials = D_T.n_one_dir_trials;
for j = 1:n_reps
    [~,~, ps, ks, ~] = D_T.basic_decode(false, [], []);
    [~,~, ps_s, ks_s, ~] = D_T.basic_decode(true, [], []);
    [~,~, ps_d, ks_d, ~] = D_T.basic_decode(false, [], [], my_algs('ecoclin', 'shuf'));
    remapper = reshape([(1:20)',(21:40)']', 1, []);
    ps = remapper(ps);
    ks = remapper(ks);
    ps_s = remapper(ps_s);
    ks_s = remapper(ks_s);
    ps_d = remapper(ps_d);
    ks_d = remapper(ks_d);
    
    C(:,:,j) = confusionmat(ks, ps); %./ (2*D_T.n_one_dir_trials);
    C_s(:,:,j) = confusionmat(ks_s, ps_s);% ./ (2*D_T.n_one_dir_trials);
    C_d(:,:,j) = confusionmat(ks_d, ps_d);
    %C_perfect = confusionmat(ks_s, ks);% ./ (2*D_T.n_one_dir_trials);
    %CDiff(:,:,j) = C_s - C;
    %CsC = (C_s - C)./C;
    %CsC(isnan(CsC)) = 0;
    progressbar(j/n_reps);
end

res.num_trials = num_trials;
%res.C = squeeze(mean(C,3));
%res.C_s = squeeze(mean(C_s,3));
res.C = C;
res.C_s = C_s;
res.C_d = C_d;
res.source = D_T.source_path;
toc
end