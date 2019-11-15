function match_list = make_matchlist(ds_list)
% A convenience function for making a 'match_list' to be used in the
% instantiation of a MultiDay object.
%
% Example usage:
%   ds_list = {12, m1d12; 13, m1d13; 14, m1d14};
%   match_list = make_matchlist(ds_list)
%
%   will perform the necessary alignment runs to yield:
%   
%   match_list = {12, 13, m_12to13, m_13to12;
%                 12, 14, m_12to14, m_14to12;
%                 13, 14, m_13to14, m_14to13};
%
%   which can then be used to instantiate a MultiDay object as:
%
%   md = MultiDay(ds_list, match_list);
%

num_ds = size(ds_list, 1);

% Allocate the match_list
num_matches = num_ds*(num_ds-1)/2;
match_list = cell(num_matches, 4);

match_idx = 1;
for i = 1:(num_ds-1)
    day_i = ds_list{i,1};
    ds_i = ds_list{i,2};
    for j = (i+1):num_ds
        day_j = ds_list{j,1};
        ds_j = ds_list{j,2};
        
        fprintf('%s: Matching day %d and day %d...\n',...
            datestr(now), day_i, day_j);
        close all;
        [m_itoj, m_jtoi] = run_alignment(ds_i, ds_j, 'fast');
        
        % Save result to match_list
        match_list{match_idx, 1} = day_i;
        match_list{match_idx, 2} = day_j;
        match_list{match_idx, 3} = m_itoj;
        match_list{match_idx, 4} = m_jtoi;
        
        match_idx = match_idx + 1;
    end
end