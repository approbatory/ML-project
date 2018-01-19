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


for i = 1:length(days)
    for j = 1:(i-1)
        mask = check_ident(days(i), days(j), matching(:,i), matching(:,j), match_list);
        matching = matching(mask,:);
    end
end
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

function res = check_ident(day_A, day_B, ix_A, ix_B, match_list)
N = length(ix_A);
if length(ix_A) ~= length(ix_B)
    error('indices vectors must have the same length');
end
table = cell2mat(match_list(:,1:2));
maps = match_list(:,3:4);
%checking for order A B
order_AB = all(table == [day_A day_B],2);
order_BA = all(table == [day_B day_A],2);
if any(order_AB)
    row = find(order_AB,1);
    col = 1;
elseif any(order_BA)
    row = find(order_BA,1);
    col = 2;
else
    error('comparison between days %d, %d not found', day_A, day_B);
end
map = maps{row,col};
res = false(size(ix_A));
for i = 1:N
    c = map{ix_A(i)};
    if isempty(c)
        res(i) = false;
    else
        res(i) = c(1) == ix_B(i);
    end
end
end