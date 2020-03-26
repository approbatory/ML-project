%parpool(20);
%parfor i = 1:107
%    res(i) = med_loadings_compute_(i);
%    fprintf('Done %d\n', i);
%end
%
%save medload_diet.mat res
%%
%mlc_ = @med_loadings_compute_;
function res = med_loadings_compute(index, n_type)
if ~exist('n_type', 'var')
    n_type = 'rawTraces';
end
o_ = DecodeTensor.cons_filt(index, true);
d = DecodeTensor(o_, n_type);
n_max = size(d.data_tensor,1);
res.n_sizes = [unique(ceil(10.^(log10(2):0.1:log10(n_max)))) n_max];
n_reps = 20;
for n_i = 1:numel(res.n_sizes)
    for i = 1:n_reps
        my_res = d.neuron_subset_noise_props(res.n_sizes(n_i));
        %if length(ml) < 50
        %    ml = [ml zeros(1, 50 - length(ml))];
        %    mls = [mls zeros(1, 50 - length(mls))];
        %else
        %    ml = ml(1:50);
        %    mls = mls(1:50);
        %end
        %res.median_loadings(i, n_i, :) = ml;
        %res.median_loadings_s(i, n_i, :) = mls;
        res.output{i, n_i} = my_res;
        %progressbar(n_i/numel(res.n_sizes), i/n_reps);
    end
end
res.mouse_name = d.mouse_name;
res.source_path = d.source_path;
%dname = 'loadings_records';
%if ~exist(dname, 'dir')
%    mkdir(dname);
%end
%
%fname = sprintf('med_loadings_%d_%s.mat', index, d.mouse_name);
%fname = fullfile(dname, fname);
%save(fname, 'n_sizes', 'median_loadings', 'median_loadings_s');
end




%{
function res = med_loadings_compute(index, n_type)
if ~exist('n_type', 'var')
    n_type = 'rawTraces';
end
o_ = DecodeTensor.cons_filt(index, true);
d = DecodeTensor(o_, n_type);
n_max = size(d.data_tensor,1);
res.n_sizes = [unique(ceil(10.^(log10(2):0.1:log10(n_max)))) n_max];
n_reps = 20;
for n_i = 1:numel(res.n_sizes)
    for i = 1:n_reps
        [ml, mls] = d.signal_loadings(res.n_sizes(n_i));
        if length(ml) < 50
            ml = [ml zeros(1, 50 - length(ml))];
            mls = [mls zeros(1, 50 - length(mls))];
        else
            ml = ml(1:50);
            mls = mls(1:50);
        end
        res.median_loadings(i, n_i, :) = ml;
        res.median_loadings_s(i, n_i, :) = mls;
        progressbar(n_i/numel(res.n_sizes), i/n_reps);
    end
end
res.mouse_name = d.mouse_name;
res.source_path = d.source_path;
%dname = 'loadings_records';
%if ~exist(dname, 'dir')
%    mkdir(dname);
%end
%
%fname = sprintf('med_loadings_%d_%s.mat', index, d.mouse_name);
%fname = fullfile(dname, fname);
%save(fname, 'n_sizes', 'median_loadings', 'median_loadings_s');
end
%}