function remainers = cell_distance_filter(distance_matrix, d_min)
remainers = true(size(distance_matrix,1),1);
target_map = (distance_matrix < d_min) & (distance_matrix ~= 0);
while true
    [i_s, j_s] = find(target_map);
    
    n_options = numel(i_s);
    chosen_option = randperm(n_options, 1);
    chosen_cell1 = i_s(chosen_option);
    chosen_cell2 = j_s(chosen_option);
    
    target_map(chosen_cell1, :) = false;
    target_map(:, chosen_cell1) = false;
    remainers(chosen_cell1) = false;
    if ~any(target_map(:))
        return;
    end
end
end

%%USE MCMC with gibbs sampling to select a UNIQUELY SIZED subset that
%%satisfies the distance requirement and has the maximal allowed number of
%%cells, i.e., minimizes Sum_i -s_i + J*Sum_j s_i*s_j*I[dist(i,j) > d & dist(i,j) ~= 0]
%%simple ising model s_i s_j A_ij = s_i A_ij s_j = s' * A * s
%%and J >> 1