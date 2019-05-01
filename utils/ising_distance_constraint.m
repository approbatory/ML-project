%%The problem of constraining the distance between cells to be greater than
%%a given value, yet maximizing the number of cells in the remaining subset
%%is identical to the minimization problem of:
%%H = -Sum_i q_i + J/2*Sum_i Sum_j q_i q_j I[dist(i,j) < d & i~=j]
%%where q has values in {0,1}. To remap to {-1,1}, define q = (s+1)/2
%%then, neglecting constants,
%%H = -1/2*Sum_i s_i + J/8*Sum_i Sum_j (s_i+1)(s_j+1)I[dist(i,j) < d & i~=j]
%%H = -1/2*Sum_i s_i + J/8*Sum_i Sum_j (s_i*s_j + s_i + s_j)I[dist(i,j) < d& i~=j]
%%easier to work with q variables
function q = ising_distance_constraint(cellAnatomicLocat, d)
%coord_x = tracesEvents.cellAnatomicLocat(:,1);
%coord_y = tracesEvents.cellAnatomicLocat(:,2);
coord_x = cellAnatomicLocat(:,1);
coord_y = cellAnatomicLocat(:,2);

distance_mat = sqrt((coord_x - coord_x').^2 + (coord_y - coord_y').^2);
n_cells = length(coord_x);

q = randi(2, n_cells, 1) - 1;

J = 2;%10;%1000;%2;
beta_max = 10;%100;%10;
%beta = 1000;
%d = 15;
interaction_mat = (distance_mat < d) & (distance_mat ~= 0);
H_func = @(q) -sum(q) + J/2*q'*interaction_mat*q;

tic;
if false
    i_times_q = interaction_mat * q;
    %%if index k flips from q_old to q_new, then
    %%dH = H_new-H_old = -(q_new-q_old) + J*(q_new-q_old)*Sum_i q_i I_ik
    %%if dH < 0, keep flip. else, keep flip with prob exp(-beta*dH)
    T = 1e5;
    for t = 1:T
        beta = t/T*beta_max;
        for k = 1:n_cells
            delta_q = 1-2*q(k);
            %delta_H = -delta_q + J*delta_q*interaction_mat(k,:)*q;
            delta_H = -delta_q + J*delta_q*i_times_q(k);
            if (delta_H < 0) || (rand < exp(-beta*delta_H))
                q(k) = ~q(k); %TODO: get rid of '*' in delta_H, use addition and subtraction of I columns, once it is changed
                i_times_q = i_times_q + (q(k)*2-1)*interaction_mat(:,k);
            end
        end
    end
end
q = mcmc_solver(J, beta_max, interaction_mat, 1e5, randi(2^16));
toc
H_final = H_func(q);
if nnz(interaction_mat(logical(q),logical(q))) == 0
    fprintf('No clashes, number of cells = %d\n', sum(q));
else
    fprintf('Clashes exist, try again\n');
end