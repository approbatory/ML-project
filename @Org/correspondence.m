function [T, vartable, rank_table] = correspondence(o, use_n50, plot_all)

load asymp_snr_values.mat
load single_cell_dp2.mat
load decoding_curves_fits.mat

[mat, sem_, varnames, mat_agg, sem_agg] = o.predictor_matrix;

if ~use_n50
    [asymp_ratio, asymp_ratio_conf] = ...
        uncertain_divide(asymp_snr, asymp_snr_conf,...
        asymp_snr_shuf, asymp_snr_shuf_conf);
    
    %[asymp_worst, asymp_worst_conf] = ...
    %    worst_divide(asymp_snr, asymp_snr_conf,...
    %    asymp_snr_shuf, asymp_snr_shuf_conf);
else
    %asymp_ratio = I0_fit * N_fit ./ (I0_fit_s * N_fit_s);
    
    if false
        [I0N, I0N_c] = uncertain_multiply(I0_fit, I0_conf, N_fit, N_conf);
        [I0N_s, I0N_c_s] = uncertain_multiply(I0_fit_s, I0_conf_s, N_fit_s, N_conf_s);
        [asymp_ratio, asymp_ratio_conf] = uncertain_divide(I0N, I0N_c, I0N_s, I0N_c_s);
    end
    
    %asymp_ratio = N_fit; asymp_ratio_conf = N_conf;
    asymp_ratio = 1./N_fit; asymp_ratio_conf = N_conf ./ N_fit.^2;
    
    if false
        chosen_size = 150;
        for i = 1:numel(sess)
            ix = find(n_sizes{i} == chosen_size, 1);
            assert(~isempty(ix));
            imse_at_size(i) = mean(imse{i}(:,ix));
            imse_at_size_conf(i) = sem(imse{i}(:,ix));
            
            imse_s_at_size(i) = mean(imse_s{i}(:,ix));
            imse_s_at_size_conf(i) = sem(imse_s{i}(:,ix));
        end
        [asymp_ratio, asymp_ratio_conf] = ...
            uncertain_divide(imse_at_size, imse_at_size_conf,...
            imse_s_at_size, imse_s_at_size_conf);
        asymp_ratio_conf = 1.96 * asymp_ratio_conf;
    end
    %[asymp_ratio, asymp_ratio_conf] = ...
    %    uncertain_divide(asymp_snr, asymp_snr_conf,...
    %    single_dp2, single_dp2_sem*1.96);
    
    %asymp_ratio = asymp_snr;
    %asymp_ratio_conf = asymp_snr_conf;
    
    %asymp_ratio = asymp_snr_shuf;
    %asymp_ratio_conf = asymp_snr_shuf_conf;
end
[asymp_ratio_agg, asymp_ratio_conf_agg] = agg(o.mouse, asymp_ratio', asymp_ratio_conf');

conf = sem_ * 1.96;
conf_agg = sem_agg * 1.96;

for j = 1:numel(varnames)
    [pearson(j,1), pearson_p(j,1)] = corr(mat(:,j), asymp_ratio(:), 'type', 'Pearson');
    [spearman(j,1), spearman_p(j,1)] = corr(mat(:,j), asymp_ratio(:), 'type', 'Spearman');
    [kendall(j,1), kendall_p(j,1)] = corr(mat(:,j), asymp_ratio(:), 'type', 'Kendall');
    [~, adjr2(j,1)] = Utils.regress_line(mat(:,j), asymp_ratio(:));
    
    [pearson_agg(j,1), pearson_p_agg(j,1)] = corr(mat_agg(:,j), asymp_ratio_agg(:), 'type', 'Pearson');
    [spearman_agg(j,1), spearman_p_agg(j,1)] = corr(mat_agg(:,j), asymp_ratio_agg(:), 'type', 'Spearman');
    [kendall_agg(j,1), kendall_p_agg(j,1)] = corr(mat_agg(:,j), asymp_ratio_agg(:), 'type', 'Kendall');
    [~, adjr2_agg(j,1)] = Utils.regress_line(mat_agg(:,j), asymp_ratio_agg(:));
end

T = table(pearson, pearson_p, spearman, spearman_p,...
    kendall, kendall_p, adjr2,...
    pearson_agg, pearson_p_agg, spearman_agg, spearman_p_agg,...
    kendall_agg, kendall_p_agg, adjr2_agg, 'RowNames', varnames);

if plot_all
    for j = 1:numel(varnames)
        figure;
        o.correlogram({mat(:,j), sem_(:,j)}, {asymp_ratio, asymp_ratio_conf./1.96}, true);
        xlabel(esc(varnames{j}));
        ylabel('Asymp ratio');
    end
end


vartable = array2table([mat, 1./asymp_ratio', asymp_ratio'], 'VariableNames', [varnames, {'N50', 'invN50'}]);
vartable = [table(mouse_names', 'VariableNames', {'Mouse'}), vartable];

table_var_names = vartable.Properties.VariableNames;
rank_table = table(mouse_names', 'VariableNames', {'Mouse'});
for i = 1:numel(table_var_names)
    name = table_var_names{i};
    if strcmp(name, 'Mouse')
        continue;
    end
    rank_table = [rank_table, table(value_rank(vartable.(name)), 'VariableNames', {name})];
end
end

function [r, rc] = agg(mouse, x, xc)
uniq_mice = unique(mouse);
n = numel(uniq_mice);
r = zeros(n, size(x,2));
rc = zeros(n, size(xc,2));
for i = 1:n
    filt = strcmp(mouse, uniq_mice{i});
    r(i,:) = mean(x(filt,:),1);
    rc(i,:) = sqrt(mean(xc(filt,:).^2,1));
end
end

function [quotient, quotient_uncertainty] = uncertain_divide(x, xc, y, yc)
quotient = x./y;
quotient_uncertainty = abs(x./y).*sqrt((xc./x).^2 + (yc./y).^2);
end

function [product, product_uncertainty] = uncertain_multiply(x, xc, y, yc)
product = x.*y;
product_uncertainty = abs(product) .* sqrt((xc./x).^2 + (yc./y).^2);
end

function [quotient, quotient_uncertainty] = worst_divide(x, xc, y, yc)
lx = x - xc;
ux = x + xc;
ly = y - yc;
uy = y + yc;

lq = lx./uy;
uq = ux./ly;

quotient = (uq + lq) ./ 2;
quotient_uncertainty = (uq - lq) ./ 2;
end