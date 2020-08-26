function res = int_time_runner(i)
[c_err, b_err] = bin_integration_time_checker(i, [2 5 10 20 30 40 50], [1 2 4 6 8 10], 20);
res.c_err = c_err;
res.b_err = b_err;
end