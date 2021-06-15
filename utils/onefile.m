function f = onefile(f)
L = filelist(f);
assert(~isempty(L), 'File %s not found', f);
assert(numel(L)==1, 'More than one file fits %s', f);
f = L{1};
end