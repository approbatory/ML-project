function res = events_nbin_checker(i_usable)

res = bin_integration_time_checker(i_usable, [2 5 10 20 30 50],...
    1, 80, my_algs('ecoclin'));