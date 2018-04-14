if exist('_data', 'dir')
    return;
end

!mkdir -p _data
!mv *.xy *.mp4 *.hdf5 _data

num_frames = length(load(file_pattern('_data', '*.xy')));
S = dir('_data/*.xy');
basename = S.name(1:end-3);
fid = fopen(fullfile('_data', [basename '.txt']), 'w');
fprintf(fid, 'east north north 1.0 1 2 3 %d', num_frames);
fclose(fid);

getfile = @(ext) fullfile('_data', file_pattern('_data', ['*.' ext], true));
maze = getfile('txt');
behavior = getfile('mp4');
tracking = getfile('xy');
miniscope = getfile('hdf5');

fid = fopen('data_sources.m', 'w');
fprintf(fid, [
'function sources = data_sources\n', ...
'\n', ...
'sources = struct(...\n',...
'\t''maze'', ''%s'',...\n',...
'\t''behavior'', ''%s'',...\n',...
'\t''tracking'', ''%s'',...\n',...
'\t''miniscope'', ''%s'');\n'
], maze, behavior, tracking, miniscope);
fclose(fid);