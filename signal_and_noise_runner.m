n_sessions = 107;
d_neu = 10;
n_reps = 20;

dm2 = cell(1, n_sessions);
sm = dm2; sms = dm2; n_sizes = dm2;
progressbar('sess', 'num neu', 'rep');
for i = 1:n_sessions
    
    d = DecodeTensor.cons_filt(i);
    [dm2{i}, sm{i}, sms{i}, n_sizes{i}] = d.signal_and_noise_descriptors_series(d_neu, n_reps);
    
    progressbar(i/n_sessions, [], []);
    
end

save signal_and_noise_save.mat dm2 sm sms n_sizes