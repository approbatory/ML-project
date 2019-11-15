function events = detect_all_events(ds)

save_to_file = true;
if ds.full_num_frames ~= ds.trial_indices(end,end)
    fprintf('Warning: Event detection should be performed for all trials, including probes!\n');
    fprintf('  Current results will not be saved to file.\n');
    save_to_file = false;
end

events = cell(ds.num_cells, 1);
for k = 1:ds.num_cells
    if (k==1) || (mod(k,50)==0)
        fprintf('%s: At cell %d of %d...\n', datestr(now), k, ds.num_cells);
    end
    events{k} = detect_events(ds, k, 'noprompt');
end

events = cell2mat(events); % Convert cell to array of structs

% Save to file
if save_to_file
    timestamp = datestr(now, 'yymmdd-HHMMSS');
    event_savename = sprintf('events_%s.mat', timestamp);
    save(event_savename, 'events', '-v7.3');
    
    ds.load_events(event_savename);
end