function [poss, err, errmap] = decode_end_svm(ds, step_size, final_pos)
MIN_POS = 0;
if ~exist('step_size', 'var')
    step_size = 0.05;
end
if ~exist('final_pos', 'var')
    final_pos = 0.4;
end

end_labels = {ds.trials.end};
evs = {ds.trials.events};
poss = MIN_POS:step_size:final_pos;
[frame_of_interest, ~] = find_frame_at_pos(ds, poss);

err = zeros(size(poss));
for j = 1:length(poss)
    %needs to be transposed since svm needs observations in rows
    X = gen_X_at_frames(evs, frame_of_interest(:,j)); 
    %SVMModel = fitcsvm(X', end_labels, 'KernelFunction', 'linear',...
    %    'Standardize', true, 'ClassNames', {'south','north'});
    CVMdl = fitclinear(X, end_labels, 'ObservationsIn', 'columns',...
        'KFold', length(end_labels), 'Learner', 'svm',...
        'ClassNames', {'south', 'north'});
    err(j) = kfoldLoss(CVMdl);
end
errmap = nan; %UNDER CONSTRUCTION
end