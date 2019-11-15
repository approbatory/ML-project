function filename = get_most_recent_file(path_to, pattern)
% Usage:
%   get_most_recent_file('ica001', 'ica_*.mat')

files = dir(fullfile(path_to, pattern));

filename = '';
if ~isempty(files)
    datenums = [files.datenum];
    [~, sorted] = sort(datenums, 'descend');

    filename = fullfile(path_to, files(sorted(1)).name);
end