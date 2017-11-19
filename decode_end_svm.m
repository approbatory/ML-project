function [poss, err, err_map] = decode_end_svm(ds, step_size, final_pos)
MIN_POS = 0;
if ~exist('step_size', 'var')
    step_size = 0.05;
end
if ~exist('final_pos', 'var')
    final_pos = 0.4;
end

label_to_class = containers.Map({'north','south','east','west'},{1,2,3,4});

end_labels = {ds.trials.end};
evs = {ds.trials.events};
poss = MIN_POS:step_size:final_pos;
[frame_of_interest, ~] = find_frame_at_pos(ds, poss);

ks = zeros(1,size(frame_of_interest,1));
K = 4;
for i=1:size(frame_of_interest,1)
    ks(i) = label_to_class(end_labels{i});
end

err = zeros(size(poss));
err_map = zeros(size(frame_of_interest));
for j = 1:length(poss)
    %needs to be transposed since svm needs observations in rows
    X = gen_X_at_frames(evs, frame_of_interest(:,j)); 
    %SVMModel = fitcsvm(X', end_labels, 'KernelFunction', 'linear',...
    %    'Standardize', true, 'ClassNames', {'south','north'});
%    CVMdl = fitclinear(X, end_labels, 'ObservationsIn', 'columns',...
%        'KFold', length(end_labels), 'Learner', 'svm',...
%        'ClassNames', {'south', 'north'});
%    err(j) = kfoldLoss(CVMdl);
    mask = false(1,size(X,2));
    
    count = 0;
    for i = 1:size(X,2) %leave one out xval
        mask(i) = true;
        
        X_train  = X(:,~mask);
        ks_train = ks(~mask);
        X_test   = X(:, mask);
        ks_test  = ks(mask);
        model = fitclinear(X_train, ks_train, 'ObservationsIn', 'columns',...
            'Learner', 'svm', 'ClassNames', [1 2]);
        ks_predicted = predict(model, X_test, 'ObservationsIn', 'columns');
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