function sherlock_aggregator_unit(source)
cd ../ML-project
addpath utils decoding
source_dir = dir(source);
saves_dir = [source '_records'];
try
    merged = Analyzer.recreate_merge(saves_dir);
    merged.save_res(source_dir.folder);
    rmdir(saves_dir, 's');
catch me
    fprintf('%s / %s\n', me.identifier, me.message);
end
exit