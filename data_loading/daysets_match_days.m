function matched_cells = daysets_match_days(daysets)
matched_cells = cell(1,numel(daysets));
for m_ix = 1:numel(daysets)
    direc = daysets{m_ix}(1).directory;
    match_file_location = fullfile(direc,'match-fix.txt');
    if ~exist(match_file_location, 'file')
        %multiday input col 1
        col1 = cellfun(@(s)str2double(s(end-1:end)), {daysets{m_ix}.day}, 'UniformOutput', false)';
        col2 = cellfun(@to_daysum, {daysets{m_ix}.directory}, {daysets{m_ix}.day}, 'UniformOutput', false)';
        ds_list = [col1 col2];
        
        load(file_pattern(fullfile(direc, '_match-fix'), 'matchlist*'), 'match_list'); %loads match_list
        md = MultiDay(ds_list, match_list);
        matched_indices = md.matched_indices;
        save(match_file_location, 'matched_indices', '-ascii');
    end
    matched_indices = load(match_file_location);
    if size(matched_indices,2) ~= numel(daysets{m_ix})
        error('cell matching does not match number of days.');
    end
    for d_ix = 1:numel(daysets{m_ix})
        %daysets{m_ix}(d_ix).res.matched_cells = matched_indices(:,d_ix);
        matched_cells{m_ix}(:,d_ix) = matched_indices(:,d_ix);
    end
end

    function ds = to_daysum(directory, name)
        start_dir = pwd;
        cd(directory);
        cd(name);
        ds = DaySummary(data_sources, 'cm01-fix');
        cd(start_dir);
    end


    function res = file_pattern(d, pat)
        S = dir(fullfile(d,pat));
        if numel(S) ~= 1
            error('There must be exactly one file matching the pattern: %s/%s', d, pat);
        end
        res = fullfile(S.folder, S.name);
    end
end