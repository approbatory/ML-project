function [X_partition, ks_partition, vals, masks] = partition_data(X_all, ks, classes)
vals = unique(classes);
X_partition = cell(1,length(vals));
ks_partition = cell(1,length(vals));
masks = cell(1,length(vals));
for i = 1:length(vals)
    masks{i} = classes == vals(i);
    X_partition{i} = X_all(masks{i},:,:);
    ks_partition{i} = ks(masks{i});
end
end