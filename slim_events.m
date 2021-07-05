sm = SessManager;
progressbar('sess...');
for i = 1:sm.num_usable
    s = sm.cons_usable(i, true);
    s = s{1};
    load(s);
    
    f = s(17:end);
    tot_file = fullfile('../slim_events', f);
    tot_dir = fileparts(tot_file);
    mkdir(tot_dir);
    
    tracesEvents = rmfield(tracesEvents, 'rawTraces');
    tracesEvents.HD_spikes = sparse(tracesEvents.HD_spikes);
    
    save(tot_file, 'tracesEvents');
    progressbar(i/sm.num_usable);
end