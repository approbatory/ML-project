function S = cpd_cluster_solutions(models, S)
% CPD_SCREE_PLOT, plots a scree plot given struct array of model fits
%
%     models = fit_cpd(data, ...)
%     CPD_SCREE_PLOT(MODELS)

if nargin == 1
    S = cpd_pairwise_similarities(models);
end

