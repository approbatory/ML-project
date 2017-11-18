function [poss, err, err_map] = decode_end_nb(ds, step_size, final_pos)
%DECODE_END_NB Summary of this function goes here
%   Detailed explanation goes here
MIN_POS = 0;
if ~exist('step_size', 'var')
    step_size = 0.05;
end
if ~exist('final_pos', 'var')
    final_pos = 0.4;
end

label_to_class = containers.Map({'north','south','east','west'},{1,2,3,4});
%class_to_label = containers.Map({1,2,3,4},{'north','south','east','west'});

end_labels = {ds.trials.end};
evs = {ds.trials.events};
poss = MIN_POS:step_size:final_pos;
[frame_of_interest, ~] = find_frame_at_pos(ds, poss);

ks = zeros(1,size(frame_of_interest,1));
K = 4;
for i=1:size(frame_of_interest,1)
    ks(i) = label_to_class(end_labels{i});
end

err = zeros(1,size(frame_of_interest,2));
err_map = zeros(size(frame_of_interest));
for j = 1:size(frame_of_interest,2)
    X = gen_X_at_frames(evs, frame_of_interest(:,j));
    mask = false(1,size(X,2));
    
    count = 0;
    for i = 1:size(X,2) %leave one out xval
        mask(i) = true;
        
        X_train  = X(:,~mask);
        ks_train = ks(~mask);
        X_test   = X(:, mask);
        ks_test  = ks(mask);
        [log_prior, log_conditional] = multinom_nb_encode(X_train, ks_train, K);
        ks_predicted = multinom_nb_decode(X_test, log_prior, log_conditional);
        mistakes = sum(ks_test ~= ks_predicted);
        %
        %if mistakes~=0
        %    fprintf('at pos %f: ', poss(j));
        %    fprintf('class %d | ', ks_test(ks_test ~= ks_predicted));
        %    disp(i);
        %end
        %
        total = length(ks_test);
        count = count + 1;
        err(j) = err(j) + (mistakes/total);
        err_map(i,j) = mistakes;
        
        mask(i) = false;
    end
    err(j) = err(j)/count;
end
end

