function load_definitions(org)

org.make_derived('lambda_normed', {'lambda'},...
    @(l)l./sum(l));

org.make_derived('corr_lambda_normed', {'corr_lambda'},...
    @(l)l./sum(l));

org.make_derived('signal_variance', {'lambda', 'loadings'},...
    @(lambda, loadings) lambda .* loadings'.^2);

org.make_derived('corr_signal_variance', {'corr_lambda', 'corr_loadings'},...
    @(lambda, loadings) lambda .* loadings'.^2);

org.make_derived('dmu_ipr', {'dmus'},...
    @(dmu) (sum(dmu.^2).^2 ./ sum(dmu.^4)));

org.make_derived('signal_density', {'dmus'},...
    @(dmu) (sum(dmu.^2).^2 ./ sum(dmu.^4)) ./ length(dmu));

org.make_derived('signal_var', {'dmus'},...
    @(dmu) var(dmu.^2));

org.make_derived('inv_signal_var', {'signal_var'},...
    @(x) 1./x);

org.make_derived('log_of_signal_var', {'signal_var'},...
    @log);

org.make_derived('signal_log_var', {'dmus'},...
    @(dmu) var(log(dmu.^2)));

org.make_derived('exp_negative_sigma2', {'signal_log_var'},...
    @(x)exp(-x));

org.make_derived('signal_var_norm', {'dmus'},...
    @(dmu) var(dmu.^2 ./ sum(dmu.^2)));

org.make_derived('signal_log_var_norm', {'dmus'},...
    @(dmu) var(log(dmu.^2 ./ sum(dmu.^2))));

org.make_derived('diag_decoder', {'dmus', 'noise_var'},...
    @(dmu, sigma2) dmu ./ sigma2);

org.make_derived('diag_decoder_density', {'diag_decoder'},...
    @(dmu) (sum(dmu.^2).^2 ./ sum(dmu.^4)) ./ length(dmu));

org.make_derived('neuron_snr', {'dmus', 'noise_var'},...
    @(dmu, sigma2) dmu ./ sqrt(sigma2));

org.make_derived('neuron_snr2', {'dmus', 'noise_var'},...
    @(dmu, sigma2) dmu.^2 ./ sigma2);

org.make_derived('snr_density', {'neuron_snr'},...
    @(dmu) (sum(dmu.^2).^2 ./ sum(dmu.^4)) ./ length(dmu));

org.make_derived('signal_perplexity_density', {'dmus'},...
    @(dmu)perplexity_density(dmu.^2));

org.make_derived('log_dmus', {'dmus'},...
    @(x)log(abs(x)));

org.make_derived('sigdens_log', {'log_dmus'},...
    @(dmu) (sum(dmu.^2).^2 ./ sum(dmu.^4)) ./ length(dmu));

org.make_derived('ipr_log', {'log_dmus'},...
    @(dmu) (sum(dmu.^2).^2 ./ sum(dmu.^4)) );

org.make_derived('cos_area', {'loadings', 'loadings_shuf'},...
    @area_between);

org.make_derived('loadings2', {'loadings'},...
    @(x)x.^2);

org.make_derived('loadings_shuf2', {'loadings_shuf'},...
    @(x)x.^2);

org.make_derived('cos2_area', {'loadings', 'loadings_shuf'},...
    @(l,ls)area_between(l.^2, ls.^2));

org.make_derived('cos2_ratio', {'loadings', 'loadings_shuf'},...
    @(l,ls) mean((l.^2) ./ (ls.^2)));

org.make_derived('cos_area_per_neuron', {'cos_area', 'num_neurons'},...
    @(a,b)a./b);

org.make_derived('cos2_area_per_neuron', {'cos2_area', 'num_neurons'},...
    @(a,b)a./b);
%%%
org.make_derived('num_evecs', {'loadings'},...
    @numel);

org.make_derived('loadings_ipr', {'loadings'},...
    @(l) sum(l.^2).^2 ./ sum(l.^4));

org.make_derived('loadings_ipr_shuf', {'loadings_shuf'},...
    @(l) sum(l.^2).^2 ./ sum(l.^4));

org.make_derived('harmonic_delta_loadings_ipr', {'loadings_ipr', 'loadings_ipr_shuf'},...
    @(a,b)1./(1./a-1./b));

org.make_derived('weighted_overlap', {'loadings', 'lambda'},...
    @(loadings, lambda) sum(lambda .* loadings'.^2) ./ sum(lambda));

org.make_derived('cos2_area_5', {'loadings'},...
    @(c)up_to(c.^2, 5));

org.make_derived('cos2_area_10', {'loadings'},...
    @(c)up_to(c.^2, 10));

org.make_derived('cos2_area_1', {'loadings'},...
    @(c)up_to(c.^2, 1));

org.make_derived('cos2_area_2', {'loadings'},...
    @(c)up_to(c.^2, 2));

org.make_derived('cos2_area_3', {'loadings'},...
    @(c)up_to(c.^2, 3));

org.make_derived('cos2_area_4', {'loadings'},...
    @(c)up_to(c.^2, 4));

org.make_derived('cos2_area_6', {'loadings'},...
    @(c)up_to(c.^2, 6));

org.make_derived('cos2_area_7', {'loadings'},...
    @(c)up_to(c.^2, 7));

org.make_derived('cos2_area_8', {'loadings'},...
    @(c)up_to(c.^2, 8));

org.make_derived('cos2_area_9', {'loadings'},...
    @(c)up_to(c.^2, 9));

org.make_derived('cos2_area_5_shuf', {'loadings_shuf'},...
    @(c)up_to(c.^2, 5));

org.make_derived('cos2_area_1_shuf', {'loadings_shuf'},...
    @(c)up_to(c.^2, 1));

org.make_derived('cos2_area_10_shuf', {'loadings_shuf'},...
    @(c)up_to(c.^2, 10));

org.make_derived('cos2_area_20', {'loadings'},...
    @(c)up_to(c.^2, 20));

org.make_derived('cos2_area_14', {'loadings'},...
    @(c)up_to(c.^2, 14));

org.make_derived('cos2_area_20_shuf', {'loadings_shuf'},...
    @(c)up_to(c.^2, 20));

org.make_derived('cos2_area_14_shuf', {'loadings_shuf'},...
    @(c)up_to(c.^2, 14));

for d = 1:40
    
    name = sprintf('delta_cos2_area_%d',d);
    unshuf = sprintf('cos2_area_%d',d);
    shuf = [unshuf '_shuf'];
    
    org.make_derived(unshuf, {'loadings'},...
    @(c)up_to(c.^2, d));
    
    org.make_derived(shuf, {'loadings_shuf'},...
    @(c)up_to(c.^2, d));
    
    org.make_derived(name, {unshuf, shuf},...
        @(a,b)a-b);
end

org.make_derived('delta_cos2_area_20', {'cos2_area_20', 'cos2_area_20_shuf'},...
    @(a,b)a-b);

org.make_derived('delta_cos2_area_14', {'cos2_area_14', 'cos2_area_14_shuf'},...
    @(a,b)a-b);

org.make_derived('delta_cos2_area_10', {'cos2_area_10', 'cos2_area_10_shuf'},...
    @(a,b)a-b);

org.make_derived('delta_cos2_area_5', {'cos2_area_5', 'cos2_area_5_shuf'},...
    @(a,b)a-b);

org.make_derived('delta_cos2_area_1', {'cos2_area_1', 'cos2_area_1_shuf'},...
    @(a,b)a-b);

org.make_derived('corr_cos2_area_10', {'corr_loadings'},...
    @(c)up_to(c.^2, 10));

org.make_derived('corr_cos2_area_10_shuf', {'corr_loadings_shuf'},...
    @(c)up_to(c.^2, 10));

org.make_derived('corr_cos_area', {'corr_loadings', 'corr_loadings_shuf'},...
    @area_between);

org.make_derived('corr_cos2_area', {'corr_loadings', 'corr_loadings_shuf'},...
    @(l,ls)area_between(l.^2, ls.^2));

org.make_derived('corr_cos2_ratio', {'corr_loadings', 'corr_loadings_shuf'},...
    @(l,ls) mean((l.^2) ./ (ls.^2)));

org.make_derived('corr_cos_area_per_neuron', {'corr_cos_area', 'num_neurons'},...
    @(a,b)a./b);

org.make_derived('corr_cos2_area_per_neuron', {'corr_cos2_area', 'num_neurons'},...
    @(a,b)a./b);
%%%
org.make_derived('intersect', {'loadings', 'loadings_shuf'},...
    @(a,b)find(a<b,1));

org.make_var_per_sess('place_field_mode_right', {'mus'},...
    @(mus) argmax(cell2mat(mus(1:20)),[],2));

org.make_var_per_sess('place_field_mode_left', {'mus'},...
    @(mus) argmax(cell2mat(mus(21:40)),[],2));

org.make_var_per_sess('mode_dist_right', {'place_field_mode_right'},...
    @(amv) abs(amv - amv'));

org.make_var_per_sess('mode_dist_left', {'place_field_mode_left'},...
    @(amv) abs(amv - amv'));

org.make_derived('noise_cov', {'evecs', 'lambda'}, ...
    @(evecs, lambda) evecs .* lambda.' * evecs.');

org.make_derived('noise_var', {'noise_cov'},...
    @diag);

org.make_derived('noise_ipr', {'noise_var'},...
    @(sigma2) (sum(sigma2).^2 ./ sum(sigma2.^2)) );

org.make_derived('noise_density', {'noise_var'},...
    @(sigma2) (sum(sigma2).^2 ./ sum(sigma2.^2)) ./ length(sigma2) );

org.make_derived('noise_perplexity_density', {'noise_var'},...
    @perplexity_density);

org.make_derived('eigen_noise_ipr', {'lambda'},...
    @(sigma2) (sum(sigma2).^2 ./ sum(sigma2.^2)) );

org.make_derived('eigen_noise_density', {'lambda'},...
    @(sigma2) (sum(sigma2).^2 ./ sum(sigma2.^2)) ./ length(sigma2) );

org.make_derived('noise_corr', {'corr_evecs', 'corr_lambda'}, ...
    @(evecs, lambda) evecs .* lambda.' * evecs.', true);
org.make_derived('noise_corr_shuf', {'corr_evecs_shuf', 'corr_lambda_shuf'},...
    @(evecs, lambda) evecs .* lambda.' * evecs.', true);

c_snr2 = 0.05;
org.make_derived('mean_coding_noise_corr', {'neuron_snr2', 'noise_corr'},...
    @(snr2, corr_mat) (sum(reshape(corr_mat(snr2>c_snr2, snr2>c_snr2),1,[]))-sum(snr2>c_snr2))./(sum(snr2>c_snr2).^2 - sum(snr2>c_snr2)));

org.make_derived('noise_cov', {'evecs', 'lambda'}, ...
    @(evecs, lambda) evecs .* lambda.' * evecs.', true);


org.make_var_per_sess('noise_corr_avg_right', {'noise_corr'},...
    @(nc) mean(cat(3,nc{1:20}), 3));

org.make_var_per_sess('noise_corr_avg_left', {'noise_corr'},...
    @(nc) mean(cat(3,nc{21:40}), 3));

org.make_var_per_sess('noise_corr_avg_right_shuf', {'noise_corr_shuf'},...
    @(nc) mean(cat(3,nc{1:20}), 3));

org.make_var_per_sess('noise_corr_avg_left_shuf', {'noise_corr_shuf'},...
    @(nc) mean(cat(3,nc{21:40}), 3));


org.make_var_per_sess('noise_cov_avg_right', {'noise_cov'},...
    @(nc) mean(cat(3,nc{1:20}), 3));

org.make_var_per_sess('noise_cov_avg_left', {'noise_cov'},...
    @(nc) mean(cat(3,nc{21:40}), 3));

org.make_var_per_sess('sig_corr_right', {'mus'},...
    @(mus) corr(cell2mat(mus(1:20)).'));

org.make_var_per_sess('sig_corr_left', {'mus'},...
    @(mus) corr(cell2mat(mus(21:40)).'));

org.make_var_per_sess('sig_cov_right', {'mus'},...
    @(mus) cov(cell2mat(mus(1:20)).'));

org.make_var_per_sess('sig_cov_left', {'mus'},...
    @(mus) cov(cell2mat(mus(21:40)).'));

org.make_var_per_sess('sig_cov_both', {'mus'},...
    @(mus) cov(cell2mat(mus).'));

org.make_var_per_sess('tot_cov_right', {'sig_cov_right', 'noise_cov_avg_right'},...
    @(sc, nc) sc + nc);

org.make_var_per_sess('tot_cov_left', {'sig_cov_left', 'noise_cov_avg_left'},...
    @(sc, nc) sc + nc);

org.make_var_per_sess('tot_corr_right', {'tot_cov_right'},...
    @(cv_) cv_ ./ sqrt(diag(cv_)*diag(cv_).'));

org.make_var_per_sess('tot_corr_left', {'tot_cov_left'},...
    @(cv_) cv_ ./ sqrt(diag(cv_)*diag(cv_).'));

org.make_var_per_sess('com_right', {'mus'},...
    @(mus) com(cell2mat(mus(1:20))));

org.make_var_per_sess('com_left', {'mus'},...
    @(mus) com(cell2mat(mus(21:40))));

org.make_var_per_sess('com_dist_right', {'com_right'},...
    @(amv) abs(amv - amv'));

org.make_var_per_sess('com_dist_left', {'com_left'},...
    @(amv) abs(amv - amv'));

org.make_derived('eigen_snr', {'lambda', 'dmus', 'evecs'},...
    @(l,dm,e) (e.'*dm).^2 ./ l);

org.make_derived('corr_eigen_snr', {'corr_lambda', 'dmus', 'corr_evecs'},...
    @(l,dm,e) (e.'*dm).^2 ./ l);

org.make_derived('full_loadings', {'dmus', 'evecs'},...
    @(dm, e) (e.'*dm).^2);

org.make_derived('corr_full_loadings', {'dmus', 'corr_evecs'},...
    @(dm, e) (e.'*dm).^2);

org.make_derived('avg_signal', {'dmus'},...
    @(dm) mean(dm.^2));

org.make_derived('total_signal', {'dmus'},...
    @(dm) sum(dm.^2));

org.make_var_per_sess('nc_by_dist_r', {'noise_corr_avg_right', 'mode_dist_right'},...
    @(nc, d)show_vars(get_pairs(d),get_pairs(nc),0,20,20));
org.make_var_per_sess('nc_by_dist_l', {'noise_corr_avg_left', 'mode_dist_left'},...
    @(nc, d)show_vars(get_pairs(d),get_pairs(nc),0,20,20));
org.make_var_per_sess('nc_by_dist_r_shuf', {'noise_corr_avg_right_shuf', 'mode_dist_right'},...
    @(nc, d)show_vars(get_pairs(d),get_pairs(nc),0,20,20));
org.make_var_per_sess('nc_by_dist_l_shuf', {'noise_corr_avg_left_shuf', 'mode_dist_left'},...
    @(nc, d)show_vars(get_pairs(d),get_pairs(nc),0,20,20));


org.make_var_per_sess('avg_nc', {'noise_corr_avg_right', 'noise_corr_avg_left'},...
    @(nc_r, nc_l) mean([get_pairs(nc_r) ; get_pairs(nc_l)]));

org.make_var_per_sess('abs_avg_nc', {'noise_corr_avg_right', 'noise_corr_avg_left'},...
    @(nc_r, nc_l) mean([abs(get_pairs(nc_r)) ; abs(get_pairs(nc_l))]));
end


function c = com(a)
    [n, b] = size(a);
    a = scale01(a);
    
    b_ind = 1:b;
    c = sum(b_ind .* a,2) ./ sum(a,2);
end


function a = scale01(a)
    [n, b] = size(a);
    for i = 1:n
        min_val = min(a(i,:));
        max_val = max(a(i,:));
        range = max_val - min_val;
        a(i,:) = a(i,:) - min_val;
        a(i,:) = a(i,:) ./ range;
    end
end


function a = argmax(varargin)
    [~, a] = max(varargin{:});
end


function r = area_between(a,b)

%int = find(a < b, 1);
int = [];

if isempty(int)
    r = sum(a-b);
else
    r = sum(a(1:int-1) - b(1:int-1));
end
end

function r = up_to(a, n)
n = min(n, numel(a));    
r = sum(a(1:n));
end


%function [edge_mids, mean_vals, sems] = show_vars(X, Y, low, high, count)
function mean_vals = show_vars(X, Y, low, high, count)
    %[~, edges] = histcounts(X);
    edges = linspace(low, high, count + 1);
    [mean_vals, sems] = avg_bins(X,Y,edges);
    edge_mids = (edges(1:end-1) + edges(2:end))/2;
end


function [mean_vals, sems] = avg_bins(X, Y, Xbins)
[X, ord] = sort(X);
Y = Y(ord);
max_X = max(X);

n_bins = numel(Xbins)-1;
[mean_vals, sems] = deal(zeros(1,n_bins));

for i = 1:n_bins
    bin_start = Xbins(i);
    bin_end = Xbins(i+1);
    
    ind_start = find(X >= bin_start, 1);
    assert(~isempty(ind_start), 'could not find start');
    if bin_end >= max_X
        ind_end = numel(X);
    else
        ind_end = find(X > bin_end, 1) - 1;
        assert(~isempty(ind_end), 'could not find end');
    end
    
    selection = Y(ind_start:ind_end);
    mean_vals(i) = mean(selection);
    sems(i) = sem(selection);
end
end