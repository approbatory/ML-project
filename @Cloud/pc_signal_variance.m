function varargout = pc_signal_variance(o)
%% calculate 1
spectrum = medify(o.lambda);
q = cell(1,40);
for i = [(1:19) (21:39)]
    q{i} = o.lambda{i}.*o.loadings{i}.^2';
end
signal_variances = medify(q);
if nargout ~= 0
    varargout{1} = spectrum;
    varargout{2} = signal_variances;
end
%% plot 1
if nargout == 0
    figure;
    loglog(spectrum, '-ob');
    hold on;
    loglog(signal_variances, '-og');
    legend 'PC variance' 'signal variance from PC'
    xlabel 'PC index'
    ylabel 'Variance (\DeltaF/F)^2'
end
%% calculate 2
med_loadings = medify(o.loadings);
med_loadings_s = medify(o.loadings_shuf);
if nargout ~= 0
    varargout{3} = med_loadings;
    varargout{4} = med_loadings_s;
end
%% plot 2
if nargout == 0
    figure;
    hold on;
    plot(med_loadings, '-ob');
    plot(med_loadings_s, '-or');
    legend 'Signal-PC loadings' 'Shuf. loadings'
    xlabel 'PC index'
    ylabel 'cos(\Delta\mu,PC_i)'
end
%% calculate 3
%because the cov matrix is singular from lack of data
%the evecs dont form a complete basis

ipr_func = @(x)sum(x.^2).^2./sum(x.^4);

signal_ipr = cellfun(ipr_func, o.loadings);
signal_ipr_s = cellfun(ipr_func, o.loadings_shuf);
spectrum_ipr = cellfun(ipr_func, o.lambda);
spectrum_ipr_s = cellfun(ipr_func, o.lambda_shuf);
if nargout ~= 0
    varargout{5} = signal_ipr;
    varargout{6} = signal_ipr_s;
    varargout{7} = spectrum_ipr;
    varargout{8} = spectrum_ipr_s;
end
%% plot 3
if nargout == 0
    figure;
    plot(spectrum_ipr, '-ob'); hold on;
    plot(spectrum_ipr_s, '-or');
    plot(signal_ipr, ':xb');
    plot(signal_ipr_s, ':xr');
    legend 'Variance IPR' 'Var. IPR shuf.' 'Signal IPR' 'Signal IPR shuf.'
    xlabel 'Bin index'
    ylabel 'IPR (# of modes)'
    set(gca, 'YScale', 'log');
end

%% calculate KL-divergence
signal_KL = cellfun(@safeKL, o.loadings, o.loadings_shuf);
N_values = cellfun(@numel, o.loadings);
%signal_perplexity = cellfun(@(n,d)n./exp(d), N_values, signal_KL);
signal_perplexity = N_values ./ exp(signal_KL);
if nargout ~= 0
    varargout{9} = signal_KL;
    varargout{10} = signal_perplexity;
end

%% calculate area between loading curves
cutoff_PC = 30;
area_cos = cellfun(@(P,Q)safe_area(P,Q,cutoff_PC), o.loadings, o.loadings_shuf);
area_cos2 = cellfun(@(P,Q)safe_area(P.^2,Q.^2,cutoff_PC), o.loadings, o.loadings_shuf);
dmu_IPR = cellfun(ipr_func, o.dmus);
num_neurons = numel(o.dmus{1});
est_sigma = @(x) std(log(x.^2));
dens_unbiased = @(x) exp(-var(log(x.^2)));
sig_dens_unbiased = cellfun(dens_unbiased, o.dmus);
log_std = cellfun(est_sigma, o.dmus);
if nargout ~= 0
    varargout{11} = area_cos;
    varargout{12} = area_cos2;
    varargout{13} = dmu_IPR;
    varargout{14} = num_neurons;
    varargout{15} = dmu_IPR ./ num_neurons;
    varargout{16} = sig_dens_unbiased;
    varargout{17} = log_std;
end
end

function a = safe_area(P, Q, cutoff)
if isempty(P) || isempty(Q)
    a = nan;
    return;
end

a = sum(P(1:min(end,cutoff)) - Q(1:min(end,cutoff)));
end

function D = safeKL(P,Q) %squares the input!
if isempty(P) || isempty(Q)
    D = nan;
    return;
end
P = P(:).'.^2; Q = Q(:).'.^2;
D = KLDiv(P,Q);
if isempty(D)
    D = nan;
end
end

function a = medify(x)
cutoff = 50;
n_bins = numel(x);
n_elems = cutoff;
a = zeros(n_bins, n_elems);
remove = cellfun(@isempty, x);
for i = 1:n_bins
    if ~remove(i)
        a(i,:) = force_size(x{i}, n_elems);
    end
end
a(remove,:) = [];
a = median(a);
end

function x = force_size(x, n)
assert(isvector(x));
[rows,~] = size(x);

if numel(x) >= n
    x = x(1:n);
else
    extra = n - numel(x);
    x = [x(:) ; zeros(extra,1)];
    if rows == 1
        x = x.';
    end
end
end