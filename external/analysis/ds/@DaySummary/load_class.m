function load_class(obj, source)

fid = fopen(source, 'r');
class = textscan(fid, '%d %s', 'Delimiter', ',');
fclose(fid);

class = class{2};
num_class = length(class);
assert(num_class == obj.num_cells,...
       sprintf('Number of labels in %s (%d) does not match number of cells (%d) in DaySummary!',...
               source, num_class, obj.num_cells));

for k = 1:obj.num_cells
    obj.cells(k).label = class{k};
end