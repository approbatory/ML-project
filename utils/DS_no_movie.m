function ds = DS_no_movie(data_sources, cm)
if isfield(data_sources, 'behavior')
    data_sources = rmfield(data_sources, 'behavior');
end
ds = DaySummary(data_sources, cm);
end