function [mat, sem, varnames, mat_agg, sem_agg] = predictor_matrix(o)

varnames = {'cos_area', 'corr_cos_area',...
    'cos2_area', 'corr_cos2_area',...
    'cos_area_per_neuron', 'corr_cos_area_per_neuron',...
    'cos2_area_per_neuron', 'corr_cos2_area_per_neuron',...
    'num_neurons', 'num_trials', ...
    'signal_density', 'sigdens_log', 'ipr_log',...
    'cos2_ratio', 'corr_cos2_ratio', 'dmu_ipr',...
    'noise_ipr', 'noise_density',...
    'eigen_noise_ipr', 'eigen_noise_density', 'avg_signal'};


[mat, sem] = deal(zeros(o.total_sessions, numel(varnames)));

for j = 1:numel(varnames)
    [mat(:,j), sem(:,j)] = o.per_sess(varnames{j});
end


uniq_mice = unique(o.mouse);
n = numel(uniq_mice);

[mat_agg, sem_agg] = deal(zeros(n, numel(varnames)));
for i = 1:n
    m = uniq_mice{i};
    filt = strcmp(o.mouse, m);
    mat_agg(i,:) = mean(mat(filt,:));
    sem_agg(i,:) = sqrt(mean(sem(filt,:).^2));
end