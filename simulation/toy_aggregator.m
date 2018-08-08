S = dir('samples');
S = S(~[S.isdir]);
num_files = numel(S);
res = cell(1,num_files);
for ix = 1:numel(S)
    s = S(ix);
    L = load(fullfile(s.folder, s.name));
    res{ix} = L.samps;
end
res = cell2mat(res);

res_mean = mean(res);
res_errb = std(res)./sqrt(length(res));

fprintf('Final result: %f +- %f\n', res_mean, res_errb);