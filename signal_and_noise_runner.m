function signal_and_noise_runner(n_, d_)
n_sessions = 107;

my_section = ceil((1:n_sessions)./n_sessions.*d_)==n_;
my_inds = find(my_section);

d_neu = 10;
n_reps = 1;%20;

dm2 = cell(1, n_sessions);
sm = dm2; sms = dm2; n_sizes = dm2;
progressbar('sess', 'num neu', 'rep');
counter = 0;
time_counter = tic;
for i = my_inds
    
    counter = counter + 1;
    d = DecodeTensor.cons_filt(i);
    session_index(counter) = i;
    [dm2{counter}, sm{counter}, sms{counter}, n_sizes{counter}] = d.signal_and_noise_descriptors_series(d_neu, n_reps);
    
    progressbar(counter/numel(my_inds), [], []);
    
end
toc(time_counter);
fprintf('Completed %d reps in this time\n', n_reps);
return;
fname = sprintf('signal_and_noise_runner_res_%d_%d.mat', n_, d_);
save(fname, 'session_index', 'dm2', 'sm', 'sms', 'n_sizes');
