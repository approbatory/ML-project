function md = create_merge_md(ds_array)
% Helper function to create MultiDay instance to represent multiple
% extraction runs on the same dataset, to be merged into one set. Used, for
% example, with iterative CNMF or creating union filter sets across days.
%
% Usage:
%   md = create_merge_md([ds1 ds2 ds3]);

num_ds = length(ds_array);
ds_list = cell(num_ds, 2);
for k = 1:num_ds
    ds = ds_array(k);
    if all(cellfun(@isempty, ds.get_class))
        fprintf('create_merge_md: DS%d is completely unlabeled. Setting all sources to be cells!\n', k);
        ds.set_labels;
    end
    
    ds_list{k,1} = k;
    ds_list{k,2} = ds;
end

match_list = cell(num_ds-1, 4);
for k = 1:(num_ds-1)
    [m_itoj, m_jtoi] = run_alignment(ds_array(k), ds_array(k+1), 'notrans', 'noprompt');
    close all;
    match_list{k,1} = k;
    match_list{k,2} = k+1;
    match_list{k,3} = m_itoj;
    match_list{k,4} = m_jtoi;
end
md = MultiDay(ds_list, match_list, 'keepall');
