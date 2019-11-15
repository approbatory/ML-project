function copy_labels(ds_source, ds_target, match)
% Copy labels from source DaySummary to target DaySummary. If no external 
% 'match' provided, copy_labels will compute the matching internally but
% assuming no translation needed

if (nargin == 2)
    match = run_alignment(ds_source, ds_target, 'notrans', 'noprompt');
end

num_labels_copied = 0;

ds_target.reset_labels();
for k = 1:ds_source.num_cells
    m = match{k};
    if ~isempty(m)
        k2 = m(1,1);
        ds_target.cells(k2).label = ds_source.cells(k).label;
        num_labels_copied = num_labels_copied + 1;
    end
end
fprintf('  %s: Copied %d labels from "%s" to "%s"\n',...
    datestr(now), num_labels_copied, inputname(1), inputname(2));