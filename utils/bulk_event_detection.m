S = dir;
top_level = pwd;
for s = S'
    disp(s.name);
    disp(s.isdir);
    if s.isdir && (s.name(1) ~= '.')
        if exist(fullfile(s.name, 'cm01-fix'), 'dir') && isempty(dir(fullfile(s.name, 'cm01-fix', 'events*.mat')))
            fprintf('Doing directory %s :\n', s.name);
            cd(s.name);
            ds = DS_no_movie(data_sources, 'cm01-fix');
            cd cm01-fix
            detect_all_events(ds);
            cd(top_level);
        end
    end
end