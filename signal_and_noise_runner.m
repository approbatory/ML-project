function res = signal_and_noise_runner(index)
%n_sessions = 107;
%n_ = index;
%d_ = n_sessions;

%my_section = ceil((1:n_sessions)./n_sessions.*d_)==n_;
%my_inds = find(my_section);

d_neu = 10;
n_reps = 1000; %6 hours

%dm2 = cell(1, n_sessions);
%sm = dm2; sms = dm2; n_sizes = dm2;
%progressbar('sess', 'num neu', 'rep');
%counter = 0;
time_counter = tic;

%counter = counter + 1;
d = DecodeTensor.cons_filt(index);
%session_index(counter) = i;
[res.dm2, res.sm, res.sms, res.n_sizes] = d.signal_and_noise_descriptors_series(d_neu, n_reps);

%progressbar(counter/numel(my_inds), [], []);

toc(time_counter);
fprintf('Completed %d reps in this time\n', n_reps);


%return;
%fname = sprintf('signal_and_noise_runner_res_%d_%d.mat', n_, d_);
%save(fname, 'session_index', 'dm2', 'sm', 'sms', 'n_sizes');
