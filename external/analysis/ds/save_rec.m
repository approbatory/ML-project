function rec_savename = save_rec(info, filters, traces) %#ok<*INUSD>

% TODO: Sanity checks on rec file structure!

% Save to file
timestamp = datestr(now, 'yymmdd-HHMMSS');
rec_savename = sprintf('rec_%s.mat', timestamp);
save(rec_savename, 'info', 'filters', 'traces', '-v7.3');
