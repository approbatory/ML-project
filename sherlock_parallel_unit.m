function sherlock_parallel_unit(random_str, source, numer, denom)
my_rand_seed = str2double(random_str);
numer = str2double(numer);
denom = str2double(denom);
rng(my_rand_seed);
fprintf('random seed %d\n\n', my_rand_seed);
addpath utils decoding
try
    %Analyzer.dispatch(source, numer, denom);
    %Analyzer.dispatch_update(source, numer, denom);
    %test_selector(str2num(source));
    DecodeTensor.dispatch_datasize(str2double(source));
catch me
    fprintf('%s / %s\n', me.identifier, me.message);
end
exit