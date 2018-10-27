%figure;
%hold on;

%c_ix = 140;

%plot(tracesEvents.rawProb(:,c_ix));
%%plot(tracesEvents.rawTraces(500:1500,100));
%%plot(tracesEvents.tresholdEvents(500:1500,100));
%plot(tracesEvents.spikeDeconvTrace(:,c_ix));
%plot(tracesEvents.spikeDeconv(:,c_ix));
%plot(tracesEvents.spikeML(:,c_ix));

%%legend rawProb rawTraces tresholdEvents spikeDeconvTrace spikeDeconv spikeML
%legend rawProb spikeDeconvTrace spikeDeconv spikeML

o = Analyzer.recreate('full_records/Analyzer_Mouse-2022-20150326_093722-linear-track-TracesAndEvents.mat_181016-144615_0.mat');
%%
my_ks = o.data.y.ks;
K = max(my_ks);
counts_y = sparse(my_ks, 1, 1);
muti_func = @(x) muti_bin(x, my_ks, K, counts_y);

muti = zeros(1,o.data.total_neurons);
for c_ix = 1:o.data.total_neurons
    sp_inds = find(o.data.X.fast_spike(:,c_ix)~=0);
    muti(c_ix) = muti_func(sp_inds);
end

tic
n_shufs = 10000;
n_neu = o.data.total_neurons;
n_samp = size(o.data.X.fast_spike,1);
muti_shuf = zeros(n_shufs, o.data.total_neurons);
n_nzspike = sum(o.data.X.fast_spike~=0);
my_ks = o.data.y.ks;
parfor sh_ix = 1:n_shufs
    for c_ix = 1:n_neu
        fake_spike = randperm(n_samp, n_nzspike(c_ix));
        muti_shuf(sh_ix, c_ix) = muti_func(fake_spike);
    end
    if mod(sh_ix,100) == 0, fprintf('%d ', sh_ix/100); end
end
toc
%%
alpha = 0.01;
pvals = mean(muti < muti_shuf);
signif_cells = pvals < alpha;

%running BH
sorted_pvals = sort(pvals);
crit = (1:numel(sorted_pvals)) ./ numel(sorted_pvals) .* alpha;
cutoff_ind = find(sorted_pvals < crit, 1, 'last');
cutoff = sorted_pvals(cutoff_ind);
bh_signif_cells = pvals <= cutoff;