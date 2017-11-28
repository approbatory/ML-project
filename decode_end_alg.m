function [poss, err, err_map, vals, masks] = decode_end_alg(preprocessor, model_generator, predictor, ds, step_size, final_pos, do_shuffle)
MIN_POS = 0;
if ~exist('step_size', 'var')
    step_size = 0.05;
end
if ~exist('final_pos', 'var')
    final_pos = 0.4;
end
if ~exist('do_shuffle', 'var')
    do_shuffle = false;
end

end_labels = {ds.trials.end};
ks = classify_labels(end_labels);

poss = MIN_POS:step_size:final_pos;
X_all = gen_all_X_at_pos_closest(ds, poss);
start_labels = {ds.trials.start};
ks_start = classify_labels(start_labels);
[X_partition, ks_partition, vals, masks] = partition_data(X_all, ks, ks_start);
err = cell(1,length(vals));
err_map = cell(1,length(vals));
for i = 1:length(vals)
    [err{i}, err_map{i}] = evaluate_on_data(X_partition{i}, ks_partition{i}, preprocessor, model_generator, predictor, do_shuffle);
end
end

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

function ks = classify_labels(labels)
label_to_class = containers.Map({'north','south','east','west'},{1,2,3,4});
ks = zeros(1, length(labels));
for i = 1:length(labels)
    ks(i) = label_to_class(labels{i});
end
end

function [err, err_map] = evaluate_on_data(X_all, ks, preprocessor, model_generator, predictor, do_shuffle)
[M,~,P] = size(X_all);
err = zeros(1,P);
err_map = zeros(M,P);

X_all = preprocessor(X_all);
for j = 1:P
    X = X_all(:,:,j);
    mask = false(1,M);
    
    count = 0;
    for i = 1:M %leave one out xval
        mask(i) = true;
        
        X_train  = X(~mask,:);
        ks_train = ks(~mask);
        if(do_shuffle) 
            X_train = shuffle(X_train, ks_train);
        end
        
        X_test   = X(mask,:);
        ks_test  = ks(mask);
        
        model = model_generator(X_train, ks_train);
        ks_predicted = predictor(model, X_test);
        mistakes = sum(ks_test ~= ks_predicted);

        total = length(ks_test);
        count = count + 1;
        err(j) = err(j) + (mistakes/total);
        err_map(i,j) = mistakes;
        
        mask(i) = false;
    end
    err(j) = err(j)/count;
end
end