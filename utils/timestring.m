function timestamp = timestring
rand_gen = rng;
seed = rand_gen.Seed;
timestamp = datestr(now, 'yymmdd-HHMMSS');
timestamp = sprintf('%s_%d', timestamp, seed);