cohort11 = auto_dayset({'c11m1', 'c11m2', 'c11m3', 'c11m5'});
matched_cells_cohort11 = daysets_match_days(cohort11);

fillings = {'traces', 'copy_zeroed', 'copy', 'binary', 'box'};

res_cohort11_matched = cell(1,numel(fillings));
res_cohort11_unmatched =  cell(1,numel(fillings));
fprintf('out of %d\n', numel(fillings));
for i = 1:numel(fillings)
    disp(i);
    res_cohort11_matched{i} = daysets_end_decoding(cohort11, ...
        'filling', fillings{i}, 'matched_cells', matched_cells_cohort11, 'mode', 'Error');
    disp('matched');
    plot_daysets_end_decoding(cohort11, res_cohort11_matched{i},...
        'save_to', ['graphs/error_decoding/cohort11/' fillings{i} '/matched'], 'suppress', true);
    res_cohort11_unmatched{i} = daysets_end_decoding(cohort11, ...
        'filling', fillings{i}, 'mode', 'Error');
    disp('unmatched');
    plot_daysets_end_decoding(cohort11, res_cohort11_unmatched{i},...
        'save_to', ['graphs/error_decoding/cohort11/' fillings{i} '/unmatched'], 'suppress', true);
end

%%
res_cohort11_matched = cell(1,numel(fillings));
res_cohort11_unmatched =  cell(1,numel(fillings));
fprintf('out of %d\n', numel(fillings));
for i = 1:numel(fillings)
    disp(i);
    res_cohort11_matched{i} = daysets_end_decoding(cohort11, ...
        'filling', fillings{i}, 'matched_cells', matched_cells_cohort11, 'mode', 'End');
    disp('matched');
    plot_daysets_end_decoding(cohort11, res_cohort11_matched{i},...
        'save_to', ['graphs/end_decoding/cohort11/' fillings{i} '/matched'], 'suppress', true);
    res_cohort11_unmatched{i} = daysets_end_decoding(cohort11, ...
        'filling', fillings{i}, 'mode', 'End');
    disp('unmatched');
    plot_daysets_end_decoding(cohort11, res_cohort11_unmatched{i},...
        'save_to', ['graphs/end_decoding/cohort11/' fillings{i} '/unmatched'], 'suppress', true);
end
