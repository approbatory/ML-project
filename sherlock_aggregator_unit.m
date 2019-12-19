function sherlock_aggregator_unit(~)
% addpath utils decoding
% source_dir = dir(source);
% saves_dir = [source '_records'];
% try
%     merged = Analyzer.recreate_merge(saves_dir);
%     merged.save_res(source_dir.folder);
%     rmdir(saves_dir, 's');
% catch me
%     fprintf('%s / %s\n', me.identifier, me.message);
% end
% exit
addpath utils decoding
try
%     table_name = 'decoding';
%     field_names = {'Mouse', 'Setting', 'NumNeurons', 'DataSize', 'MeanErrors', 'MSE'};
%     S = dir('records');
%     dbfile = 'decoding.db';
%     if ~exist(dbfile, 'file')
%         conn = sqlite(dbfile, 'create');
%         conn.exec('CREATE TABLE decoding(Mouse text, Setting text, NumNeurons int, DataSize int, MeanErrors real, MSE real);');
%     else
%         conn = sqlite(dbfile); 
%     end
%     for i = 1:numel(S)
%         if ~S(i).isdir
%             load(fullfile(S(i).folder, S(i).name));
%             for j = 1:numel(db_queue)
%                 conn.insert(table_name, field_names, db_queue{j});
%             end
%         end
%     end
%     conn.close;
    DecodeTensor.aggregate_results('db_file', 'decoding_all_sess_IED.db');
catch me
    fprintf('%s / %s\n', me.identifier, me.message);
end
exit
