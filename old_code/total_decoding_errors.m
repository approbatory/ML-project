%Pablo asked about the raw decoding score, for each mouse:

%using bins that are 5.9cm wide and an error comparing true bin center vs
%predicted bin center: we have MSE and mean error:

dbfile = 'decoding_with_mindist.db';

conn = sqlite(dbfile); %remember to close it

num_mice = 10;
neu_size = cell(2,num_mice);
errs = cell(1, num_mice);

fprintf('\tMouse\t\tSession\tNumNeurons\tDataSize\tMeanError\tMSE\n');
for i = 1:num_mice
    [~, mouse_name, session_id] = DecodeTensor.default_datasets(i);
    
    neu_size(:,i) = conn.fetch(...
        ['select max(NumNeurons), max(DataSize) from decoding where SessionID = ''' session_id ''';']);
    errs{i} = conn.fetch(...
        ['select MeanErrors, MSE from decoding where SessionID = ''' session_id ''' and NumNeurons = ' num2str(neu_size{1,i}) ...
        ' and DataSize = ' num2str(neu_size{2,i}) ' and MinDist = 0.0;']);
    
    m_errs = mean(cell2mat(errs{i}));
    e_errs = std(cell2mat(errs{i})) ./ sqrt(size(errs{i},1));
    
    fprintf('\t%s\t%s\t%d\t\t%d\t\t%.3f+-%.3f\t%.3f+-%.3f\n', mouse_name, session_id,...
        neu_size{1, i}, neu_size{2, i}, m_errs(1), e_errs(1), m_errs(2), e_errs(2));
end
conn.close;