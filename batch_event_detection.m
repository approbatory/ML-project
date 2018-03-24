super = '..';
directory = 'c11m1';
S = dir(fullfile(super, directory, [directory '*']));

for s = S(5:7)'
    cd(s.folder);
    cd(s.name);
    if ~exist('data_sources.m', 'file')
        fid = fopen('data_sources.m', 'w');
        fprintf(fid, 'function sources = data_sources()\n\nsources = struct(...\n\t''maze'',\t''_data/%s_ti2.txt'',...\n\t''behavior'',  ''_data/%s_ti2.mp4'',...\n\t''tracking'',  ''_data/%s_ti2.xy'');\n', s.name, s.name, s.name);
        fclose(fid);
    end
    ds = DaySummary(data_sources, 'cm01-fix');
    cd cm01-fix
    fprintf('Commence event detection for %s: ', s.name);
    detect_all_events(ds);
    fprintf(' .... Done\n');
end