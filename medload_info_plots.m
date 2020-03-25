function medload_info_plots(res)
% The variable res needs to be defined. To get it run
% med_loadings_compute_(70) or any other index for a session
p = Pub(21, 4*2, 'rows', 2, 'columns', 3);


n_bins = 20;
n_dirs = 2;
n_comps = 50;
%full cells
out = res.output{1,end};

bins = meshgrid(1:n_bins,1:n_dirs);
bins = arrayfun(@(x) x*ones(n_comps,1), bins, 'UniformOutput', false);

eigvals = out.noise_spectrum(:, 2:n_bins);
loadings = out.el_pre(:, 2:n_bins);
bins = bins(:, 2:n_bins);

eigvals = cellfun(@(x)x(1:50), eigvals, 'UniformOutput', false);
loadings = cellfun(@(x)x(1:50)', loadings, 'UniformOutput', false);
indices = cellfun(@(x)(1:numel(x))', eigvals, 'UniformOutput', false);

%%
use_all = true;
specific_dir = 2;
specific_bin = 11;
if use_all
    eigvals = cat(1, eigvals{:});
    loadings = cat(1, loadings{:});
    indices = cat(1, indices{:});
    bins = cat(1, bins{:});
else
    eigvals = eigvals{specific_dir, specific_bin};
    loadings = loadings{specific_dir, specific_bin};
    indices = indices{specific_dir, specific_bin};
    bins = bins{specific_dir, specific_bin}
end
%%

%figure(1);
p.panel(1, 'xlab', 'Correlation eigenvalue', 'ylab', 'max_i|cos(PC_i, \Delta\mu)|', 'title', sprintf('%s, one session, all bins', res.mouse_name));
up_to = 7;
scatter(eigvals(indices > up_to), abs(loadings(indices > up_to)), 1, 'k');
hold on;
scatter(eigvals(indices <= up_to), abs(loadings(indices <= up_to)), 1, indices(indices <= up_to));
xlabel 'Correlation eigenvalue'
ylabel 'max_i|cos(PC_i, \Delta\mu)|'

set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
colormap(gca, lines(up_to));
h = colorbar;
h.Label.String = 'Fluctuation mode';
title(sprintf(' %s, one session, all bins', res.mouse_name));
%figure_format('boxsize', [3 2]/2);

%{
subplot(1,2,2);
scatter(eigvals, abs(loadings), 4, bins);
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
colormap(gca, 'hsv');
colorbar
%}
%%
up_to = 50;
cutoff_n = 100;
[n, s_mean, s_sem] = get_eigenvalue_series(res, up_to);
n_cell = cellfun(@(x,n) 0*x+n, s_mean, num2cell(n), 'UniformOutput', false);

n_vec = cat(1,n_cell{:});
s_mean_vec = cat(1,s_mean{:});
s_sem_vec = cat(1,s_sem{:});

%figure;
%errorbar(n_vec, s_mean_vec, s_sem_vec, 'o');

%figure(2);
p.panel(2, 'xlab', 'Number of cells', 'ylab', 'Correlation eigenvalue');
%figure(3);
%set(gca, 'ColorOrder', parula(50));
my_colors = jet(up_to);
curve_set = cell(1,up_to);
for c_i = 1:up_to
    for n_i = 1:numel(n)
        num_cells = n(n_i);
        s = s_mean{n_i};
        ss = s_sem{n_i};
        num_avail = numel(s);
        if c_i > num_avail
            continue;
        end
        curve_set{c_i} = [curve_set{c_i} ; num_cells s(c_i) ss(c_i)];
    end
    
    n_cells = curve_set{c_i}(:, 1)';
    n_cells = curve_set{c_i}(n_cells >= cutoff_n, 1)';
    s_val = curve_set{c_i}(n_cells >= cutoff_n, 2)';
    [slope(c_i), slope_conf(c_i)] = Utils.fitaline(n_cells, s_val);
    [expo(c_i), expo_conf(c_i)] = Utils.fitaline(log(n_cells), log(s_val));
    %figure(2);
    errorbar(curve_set{c_i}(:,1), curve_set{c_i}(:,2), curve_set{c_i}(:,3), 'Color', my_colors(c_i,:));
    hold on;
    %figure(3);
    %%errorbar(curve_set{c_i}(:,1), curve_set{c_i}(:,2) / slope(c_i), curve_set{c_i}(:,3) / slope(c_i), 'Color', my_colors(c_i,:));
    %plot(curve_set{c_i}(:,1), curve_set{c_i}(:,2) / curve_set{c_i}(end,2), 'Color', my_colors(c_i,:));
    hold on;
end
for f_i = 2:2
    %figure(f_i);
    xlabel 'Number of cells'
    if f_i == 2
        ylabel 'Correlation eigenvalue'
    else
        ylabel(sprintf('Correlation eigenvalue\n(as fraction of max)'));
    end
    %set(gca, 'XScale', 'log');
    %set(gca, 'YScale', 'log');
    %grid on;
    colormap(gca, 'jet');
    h = colorbar;
    
    h.Label.String = 'Fluctuation mode';
    %axis equal;
    %axis square;
    %figure_format('boxsize', [3 2]/2);
    h.TickLabels = (0:0.2:1)*50;
    xlim([0 500]);
end



p.panel(3, 'xlab', 'Number of cells', 'ylab', sprintf('Correlation eigenvalue\n(as fraction of max)'));
%figure(3);
%set(gca, 'ColorOrder', parula(50));
my_colors = jet(up_to);
curve_set = cell(1,up_to);
for c_i = 1:up_to
    for n_i = 1:numel(n)
        num_cells = n(n_i);
        s = s_mean{n_i};
        ss = s_sem{n_i};
        num_avail = numel(s);
        if c_i > num_avail
            continue;
        end
        curve_set{c_i} = [curve_set{c_i} ; num_cells s(c_i) ss(c_i)];
    end
    
    n_cells = curve_set{c_i}(:, 1)';
    n_cells = curve_set{c_i}(n_cells >= cutoff_n, 1)';
    s_val = curve_set{c_i}(n_cells >= cutoff_n, 2)';
    [slope(c_i), slope_conf(c_i)] = Utils.fitaline(n_cells, s_val);
    [expo(c_i), expo_conf(c_i)] = Utils.fitaline(log(n_cells), log(s_val));
    %figure(2);
    %errorbar(curve_set{c_i}(:,1), curve_set{c_i}(:,2), curve_set{c_i}(:,3), 'Color', my_colors(c_i,:));
    %hold on;
    %figure(3);
    %errorbar(curve_set{c_i}(:,1), curve_set{c_i}(:,2) / slope(c_i), curve_set{c_i}(:,3) / slope(c_i), 'Color', my_colors(c_i,:));
    plot(curve_set{c_i}(:,1), curve_set{c_i}(:,2) / curve_set{c_i}(end,2), 'Color', my_colors(c_i,:));
    hold on;
end
for f_i = 3:3
    %figure(f_i);
    xlabel 'Number of cells'
    if f_i == 2
        ylabel 'Correlation eigenvalue'
    else
        ylabel(sprintf('Correlation eigenvalue\n(as fraction of max)'));
    end
    %set(gca, 'XScale', 'log');
    %set(gca, 'YScale', 'log');
    %grid on;
    colormap(gca, 'jet');
    h = colorbar;
    
    h.Label.String = 'Fluctuation mode';
    %axis equal;
    %axis square;
    %figure_format('boxsize', [3 2]/2);
    h.TickLabels = (0:0.2:1)*50;
    xlim([0 500]);
end


%{
figure;
errorbar(1:up_to, expo, expo_conf);
xlabel 'Fluctuation mode'
ylabel 'Eigenvalue growth exponent'
%}
p.panel(4, 'xlab', 'Bin of signal direction',...
    'ylab', 'Bin of noise cov. matrix',...
    'title', 'Noise along signal');
C = cross_noise(res, 'noise_cov');
C = C(:,2:end);
imagesc(C);
axis image;
h = colorbar;
h.Label.String = '\Delta\mu^TC\Delta\mu / ||\Delta\mu||^2';


p.panel(5, 'xlab', 'Bin of signal direction',...
    'ylab', 'Bin of noise corr. matrix',...
    'title', 'Normalized noise along signal');
C = cross_noise(res, 'noise_corr');
C = C(:,2:end);
imagesc(C);
axis image;
h = colorbar;
h.Label.String = '\Delta\mu^TC\Delta\mu / ||\Delta\mu||^2';


p.format;
p.print('supplements_pdf', sprintf('eigenvalue_growth_%s', res.mouse_name));
end
%%
%ges_ = @get_eigenvalue_series;
function [n,s_mean,s_sem] = get_eigenvalue_series(res, max_val)
n = res.n_sizes;
out = res.output;
%out = cellfun((@(x)cellfun(@(y)y(1:min(max_val, numel(y))), x.noise_spectrum, 'UniformOutput', false)), out, 'UniformOutput', false);
out = Utils.cf_(@(x)Utils.cf_(@(y)y(1:min(max_val,numel(y))),x.noise_spectrum),out);
for n_i = 1:numel(n)
    %num_neurons = n(n_i);
    x = out(:, n_i);
    
    x = cellfun(@(y) median([y{:}],2), x, 'UniformOutput', false);
    x = [x{:}];
    s_mean{n_i} = mean(x,2);
    s_sem{n_i} = std(x,0,2) ./ sqrt(size(x,2));
end
end

function C = cross_noise(res, fieldname)
b = 20;
samp = size(res.output,1);
C = zeros(b, b, samp);
out = res.output(:,end);
for s_i = 1:samp
    corr_mats = out{s_i}.(fieldname);
    for b_i = 1:b
        corr_mat = corr_mats{1, b_i};
        for b_j = 2:b
            sig_direc = out{s_i}.sig_direc_pre{1, b_j};
            C(b_i, b_j, s_i) = sig_direc.' * corr_mat * sig_direc;
            assert(abs(norm(sig_direc) - 1) < 1e-4);
        end
    end
end
C = mean(C,3);
end