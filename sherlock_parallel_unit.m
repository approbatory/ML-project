function sherlock_parallel_unit(random_str, source)%, is_padded, distance_cutoff)
my_rand_seed = str2double(random_str);
%is_padded = str2double(is_padded);
%distance_cutoff = str2double(distance_cutoff);
rng(my_rand_seed);
fprintf('random seed %d\n\n', my_rand_seed);
setup
DecodeTensor.linear_track_path('/oak/stanford/groups/mschnitz/hazon/lintrack_home');
try
    %Analyzer.dispatch(source, numer, denom);
    %Analyzer.dispatch_update(source, numer, denom);
    %test_selector(str2num(source));
    %DecodeTensor.dispatch_datasize(str2double(source));
    %DecodeTensor.dispatch(str2double(source), is_padded, distance_cutoff);
    %DecodeTensor.dispatch_filt(str2double(source), 'IED');
    %%DecodeTensor.dispatch_datasize_filt(str2double(source));
    %SessManager.dispatch_base(str2double(source));
    sm = SessManager;
    sm.dispatch_usable(str2double(source), 'HD');
    %running only 30 trials
    %sm.dispatch_usable_trial_restricted(str2double(source));
catch me
    fprintf('%s / %s\n', me.identifier, me.message);
end
exit