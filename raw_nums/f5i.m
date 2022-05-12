function f5i
%% Figure 5i
% Depends on files: default_store.mat
% Raw numbers for:
% - ScatterX, ScatterY, ErrorBarsX, ErrorBarsY, Mouse
rng(0);

data = corr_nsv_asnr_ratio;
make_xlsx(data, 'f5i');
end

%% helper functions
function data = corr_nsv_asnr_ratio
org = Org;
org.init;

figure;
data = my_correlogram(org, 'signal_var_norm', 'asymp_ratio', true, true, true);
axis square
set(gca, 'XTick', [1 2 3]*1e-4);
Utils.fix_exponent(gca, 'x', 0);
ylabel 'Real asymp. SNR / Shuffled asymp. SNR'
xlabel 'Normalized signal variance'
end

function [data, metrics] = my_correlogram(o, var1, var2, aggregate, highqual, conf_filt)
if ~exist('conf_filt', 'var')
    conf_filt = false;
end
%if ~exist('low_bounds', 'var')
%    low_bounds = Org.default_low_bounds;
%end
v1 = o.sess_prop.(var1);
v1_conf = o.sess_prop_conf.(var1);

v2 = o.sess_prop.(var2);
v2_conf = o.sess_prop_conf.(var2);

my_mice = o.mouse;

%filt = o.sess_prop.N50 > 2*o.sess_prop_conf.N50;
%filt = filt & (v1 > 2*v1_conf) & (v2 > 2*v2_conf); %<>
%filt = o.sess_prop.num_neurons >= low_bounds(1) &...
%    o.sess_prop.num_trials >= low_bounds(2);
if ~exist('highqual', 'var')
    highqual = false;
end

if highqual
    filt = SessManager.highqual_filt_from_usable;
else
    filt = true(1,o.total_sessions);
end
if conf_filt
    filt = filt & (v1' > v1_conf') & (v2' > v2_conf');
end
if strcmp(var2, 'asymp_ratio')
    filt = filt & (v2' < 7);
    fprintf('Special for asymp_ratio, removing outliers > 7\n');
end
if any(~filt)
    fprintf('Using only %d out of %d sessions\n',...
        sum(filt), numel(filt));
    %excluded = find(~filt);
    %fprintf('Throwing out these sess indices:\n');
    %disp(excluded);
end
v1 = v1(filt);
v1_conf = v1_conf(filt);
v2 = v2(filt);
v2_conf = v2_conf(filt);
my_mice = my_mice(filt);

if aggregate
    func = @PanelGenerator.plot_regress_averaged;
    [v1_agg, v1_conf_agg] = Org.agg(my_mice, v1, v1_conf);
    [v2_agg, v2_conf_agg] = Org.agg(my_mice, v2, v2_conf);
    [p,pp,s,sp,k,kp,adjr2] = Org.corr_check(v1_agg, v2_agg);
    fprintf('Mouse-aggregated correlations %s vs. %s: adj. R^2 = %.3f\n', var1, var2, adjr2);
    fprintf('Pearson: %.3f, p = %e, %s\n', p, pp, Utils.pstar(pp));
    fprintf('Spearman: %.3f, p = %e, %s\n', s, sp, Utils.pstar(sp));
    fprintf('Kendall: %.3f, p = %e, %s\n', k, kp, Utils.pstar(kp));
    
    data.ScatterX = v1_agg;
    data.ScatterY = v2_agg;
    data.ErrorBarsX = v1_conf_agg;
    data.ErrorBarsY = v2_conf_agg;
    data.Mouse = unique(my_mice);
else
    func = @PanelGenerator.plot_regress;
    [p,pp,s,sp,k,kp,adjr2] = Org.corr_check(v1, v2);
    fprintf('Sessionwise correlations %s vs. %s: adj. R^2 = %.3f\n', var1, var2, adjr2);
    fprintf('Pearson: %.3f, p = %e, %s\n', p, pp, Utils.pstar(pp));
    fprintf('Spearman: %.3f, p = %e, %s\n', s, sp, Utils.pstar(sp));
    fprintf('Kendall: %.3f, p = %e, %s\n', k, kp, Utils.pstar(kp));
    
    data.ScatterX = v1(:);
    data.ScatterY = v2(:);
    data.ErrorBarsX = v1_conf(:);
    data.ErrorBarsY = v2_conf(:);
    data.Mouse = my_mice;
end
adjr2_ = func(v1(:), v2(:), v1_conf(:), v2_conf(:), my_mice, 'k', 'dotsize', 10);

assert(abs(adjr2 - adjr2_) < 0.01, 'Mismatch in adj. R^2 values');
xlabel(esc(var1));
ylabel(esc(var2));

metrics.p = p; metrics.pp = pp; metrics.s = s; metrics.sp = sp;
metrics.k = k; metrics.kp = kp; metrics.adjr2 = adjr2;
metrics.aggregated = aggregate;
end