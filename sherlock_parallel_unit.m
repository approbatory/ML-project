function sherlock_parallel_unit(random_str, source, numer, denom)
my_rand_seed = str2num(random_str);
numer = str2num(numer);
denom = str2num(denom);
rng(my_rand_seed);
fprintf('random seed %d\n\n', my_rand_seed);
cd ../ML-project
addpath utils decoding
try
    Analyzer.dispatch(source, numer, denom);
catch me
    fprintf('%s / %s\n', me.identifier, me.message);
end
exit