function [res, dim_names, key] = eval_one_sess(setting, func, dt)

n_neurons = setting.num_neurons_list(dt);
n_trials = setting.num_trials_list(dt);

res = cell(numel(n_neurons), numel(n_trials), setting.trial_split_reps, 2*dt.n_bins);

for n_i = 1:numel(n_neurons)
    N = n_neurons(n_i);
    for t_i = 1:numel(n_trials)
        T = n_trials(t_i);
        trials_per_part = floor(T / setting.trial_split);
        %remainder = mod(T, setting.trial_split);
        
        partition = zeros(1,dt.n_min_trials);
        for j = 1:setting.trial_split
            partition((j-1)*trials_per_part + 1: j*trials_per_part) = j;
        end
        
        if strcmp(setting.input_mode, 'tn')
            for b_i = 1:2*dt.n_bins
                if setting.adjacent_bins && (b_i == dt.n_bins || b_i == 2*dt.n_bins)
                    continue;
                end
                for r_i = 1:setting.trial_split_reps
                    partition = partition(randperm(length(partition)));
                    
                    neuron_subset = randperm(dt.total_neurons) <= N;
                    X_nt = dt.get_bin_resp(b_i);
                    X_nt = X_nt(neuron_subset,:);
                    if setting.adjacent_bins
                        X_nt_adj = dt.get_bin_resp(b_i+1);
                        X_nt_adj = X_nt_adj(neuron_subset,:);
                    end
                    
                    X_part = cell(1, setting.trial_split);
                    for j = 1:setting.trial_split
                        if setting.adjacent_bins
                            X_part{j} = {X_nt(:,partition==j).', X_nt_adj(:,partition==j).'};
                        else
                            X_part{j} = X_nt(:,partition==j).';
                        end
                    end
                    func_output = func(X_part{:});
                    res{n_i, t_i, r_i, b_i} = func_output;
                end
            end
        else
            error('unimplemented, TODO');
        end
    end
end

key.neurons = n_neurons;
key.trials = n_trials;
key.reps = setting.trial_split_reps;
key.bins = 2*dt.n_bins;

dim_names = {'neurons', 'trials', 'reps', 'bins'};
dim_sizes = size(res);

res = squeeze(res);
dim_names = dim_names(dim_sizes > 1);
end