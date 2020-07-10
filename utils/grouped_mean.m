function [m, err, varargout] = grouped_mean(x, g)

%compute m and err according to ML and guess the form for err
%my guess is based on var(means_count)/#G, maybe with a contribution
% of var(deviations)/N
%TODO::

%if ~use_model
    
    G = unique(g);
    for i = 1:numel(G)
        x_groups{i} = x(g==G(i));
    end
    
    group_means = cellfun(@mean, x_groups);
    m = mean(group_means);
    %m = mean(x);
    
    group_deviations = cellfun(@(x)x-mean(x), x_groups,...
        'UniformOutput', false);
    deviations = cell2mat(group_deviations(:));
    
    mean_counts = cellfun(@(x)x.*0+mean(x), x_groups,...
        'UniformOutput', false);
    mean_counts = cell2mat(mean_counts(:));
    
    between_var = var(mean_counts)/numel(G);
    within_var = var(deviations)/numel(x);
    
    err = sqrt(between_var + within_var);
    
%else
    ds = dataset(x,g);
    mdl = fitlme(ds, 'x ~ 1 + (1|g)');
    m_mdl = mdl.Coefficients.Estimate;
    err_mdl = mdl.Coefficients.SE;
    varargout{1} = mdl;
    
    fprintf('m abs diff:\t%e\n', abs(m - m_mdl)./m_mdl);
    fprintf('err abs diff:\t%e\n', abs(err - err_mdl)./err_mdl);
end