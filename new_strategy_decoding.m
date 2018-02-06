[dayset, matching] = my_daysets('c11m1');
for i = 1:numel(dayset)
    %[ds(i), X(i), ks(i), ~] = load_day(dayset(i), 'ds', {'deprobe', 'nocells'}, 'data', {'filling', 'binary', 
    ds(i) = quick_ds(fullfile(dayset(i).directory, dayset(i).day), 'deprobe', 'nocells');
    if i == 1
        t = zeros(1,ds(i).num_trials);
    elseif i == 3
        t = ones(1,ds(i).num_trials);
    elseif i == 2
        t = [zeros(1,50) ones(1,ds(i).num_trials-50)];
    end
    [X{i}, ks{i}, ~] = ds_dataset(ds(i), 'filling', 'binary', 'selection', 0.3,...
        'trials', strcmp({ds(i).trials.start}, 'east') & strcmp({ds(i).trials.end}, 'north'),...
        'target', num2cell(t));
end

%TODO run the decoding with the data