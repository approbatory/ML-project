function res = file_pattern(d, pat, justname)
if ~exist('justname', 'var')
    justname = false;
end
S = dir(fullfile(d,pat));
if numel(S) ~= 1
    error('There must be exactly one file matching the pattern: %s/%s', d, pat);
end
if justname
    res = S.name;
else
res = fullfile(S.folder, S.name);
end
end