function matching = find_full_match(match_list)
day_labels = cell2mat(match_list(:,1:2));
mappings = match_list(:,3:4);
days = unique(reshape(day_labels,1,[]));
ind = find(day_labels == days(1), 1);
N = length(mappings{ind});
matching = zeros(N, length(days));
matching(:,1) = (1:N)';
for i = 1:(length(days)-1)
    fst = days(i);
    snd = days(i+1);
    comp_ind = find(all(cell2mat(match_list(:,1:2)) == [fst snd],2),1);
    new_map = match_list{comp_ind,3};
    last_col = matching(:,i);
    matching(:,i+1) = fill_map(last_col, new_map);
end
matching = matching(all(matching~=0,2),:);
end


function res = fill_map(last_col, new_map)
res = zeros(size(last_col));
for ind = 1:length(last_col)
    i = last_col(ind);
    if i == 0
        continue;
    end
    c = new_map{i};
    if ~isempty(c)
        res(ind) = c(1);
    end
end
end