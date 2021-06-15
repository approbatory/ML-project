function L = filelist(L)
    if ~iscell(L) && ischar(L)
        L = glob(L);
    else
        error('Invalid input. Enter a cell or a pattern.');
    end
end