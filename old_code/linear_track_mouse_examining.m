S = dir;
off = 0;
recs = cell(7, 0);
durations = [];
neuronses = [];
for i_s = 1:numel(S)
    s = S(i_s);
    if s.name(1) == '.'
        continue;
    end
    S_sub = dir(s.name);
    for i_sub = 1:numel(S_sub)
        sub = S_sub(i_sub);
        if sub.name(1) == '.'
            continue;
        end
        if strcmp(sub.name, 'Mouse-2029-20150301_165938-linear-track-TracesAndEvents.mat')
            continue;
        end
        recs = [recs, split(sub.name, {'-', '_'})];
        load(fullfile(s.name,sub.name));
        [duration, neurons] = size(tracesEvents.rawTraces);
        durations = [durations ; duration];
        neuronses = [neuronses ; neurons];
        off = off + 1;
        fprintf('%d ', off);
        if true
            figure;
            plot(tracesEvents.position(:,1));
            title(sub.name);
            pause;
        end
    end
    %recs(i_s + off,:) =  split(s.name, '-');
end
fprintf('\n');
%%
row = 4;
for i = 1:size(recs,2)
    fprintf('''%s\n', recs{row, i});
end

%%
for i = 1:numel(durations)
    fprintf('%d\n', durations(i));
end
%%
for i = 1:numel(neuronses)
    fprintf('%d\n', neuronses(i));
end

%added mice are: 2011, 2021, 2025, 2029