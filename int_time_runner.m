function res = int_time_runner(i)
%[c_err, b_err] = bin_integration_time_checker(i, [2 5 10 20 30 40 50 60], [1 2 4 6 8 10], 20);
%res.c_err = c_err;
%res.b_err = b_err;

cluster = get_cluster('mschnitz,hns,normal', '16G', '2:00:00', 4);

n_bins_arr = [2 5 10 20 30 40 50 60];
frames_int_arr = [1 2 4 6 8 10];

n_reps = 20;

jobs = cell(numel(n_bins_arr), numel(frames_int_arr), n_reps);
for nb_i = 1:numel(n_bins_arr)
    for fi_i = 1:numel(frames_int_arr)
        for r_i = 1:n_reps
            jobs{nb_i, fi_i, r_i} = batch(cluster, @bin_integration_time_checker, 2, {i,...
                n_bins_arr(nb_i), frames_int_arr(fi_i), 1});
        end
    end
end



for nb_i = 1:numel(n_bins_arr)
    for fi_i = 1:numel(frames_int_arr)
        for r_i = 1:n_reps
            wait(jobs{nb_i, fi_i, r_i});
            outputs = fetchOutputs(jobs{nb_i, fi_i, r_i});
            res.c_err(nb_i, fi_i, r_i) = outputs{1};
            res.b_err(nb_i, fi_i, r_i) = outputs{2};
        end
    end
end
end