function ks = classify_labels(labels)
label_to_class = containers.Map({'north','south','east','west'},{1,2,3,4});
ks = zeros(1, length(labels));
for i = 1:length(labels)
    ks(i) = label_to_class(labels{i});
end
end