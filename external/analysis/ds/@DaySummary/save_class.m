function save_class(obj, outfile)

if ~exist('outfile', 'var')
    outfile = sprintf('class_%s.txt', datestr(now, 'yymmdd-HHMMSS'));
end

num_candidates = obj.num_cells;

fid = fopen(outfile, 'w');
for k = 1:num_candidates
    fprintf(fid, '%d, %s\n', k, obj.cells(k).label);
end
fclose(fid);
fprintf('%s: Class file saved to "%s"\n', datestr(now), outfile);