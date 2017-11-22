function [poss, err, err_map] = decode_end_alg(preprocessor, model_generator, predictor, ds, step_size, final_pos, do_shuffle)
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


label_to_class = containers.Map({'north','south','east','west'},{1,2,3,4});
end_labels = {ds.trials.end};
poss = MIN_POS:step_size:final_pos;
ks = zeros(1, ds.num_trials);

for i=1:ds.num_trials
    ks(i) = label_to_class(end_labels{i});
end

err = zeros(size(poss));
err_map = zeros(ds.num_trials, length(poss));

X_all = gen_all_X_at_pos_closest(ds, poss);
X_all = preprocessor(X_all);
for j = 1:length(poss)
    X = X_all(:,:,j); 
    mask = false(1,ds.num_trials);
    
    count = 0;
    for i = 1:ds.num_trials %leave one out xval
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