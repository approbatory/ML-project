function res = spike_detect(index)
dt = 1/20;
path_info = DecodeTensor.cons_filt(70);
load(path_info{1});
res.spikes = mlspike_yescal(tracesEvents.rawTraces(:, index), dt);

end