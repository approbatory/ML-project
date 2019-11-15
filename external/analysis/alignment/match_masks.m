function [match_1to2, match_2to1, M] = match_masks(masks1, masks2, ds1, ds2, varargin)
% Computes the overlaps between two sets of logical masks. Output format is
% as follows:
%
% match_1to2{k} is a [Nx2] matrix where N is the number of masks from 
%   source 2 that match mask k of source 1. The first column of the matrix 
%   is the index of the matching mask in source 2; the second column is
%   second column the overlap measure between the two masks. The list is
%   sorted in order of descending overlap measure.
%
% Note that the 'masksX' inputs may have been spatially transformed. For
% this reason, we don't simply pull the 'masks' from the corresponding
% input DaySummary objects.

% "Fast matching" mode uses two shortcuts. Given a mask i from masks1,
% if a mask is found in masks2 that exceeds a certain threshold, then:
% 1) Don't perform exhaustive search of the remaining masks2;
% 2) The matched mask in masks2 is now unavailable for matching
%    for against remaining masks1.
use_fast_matching = 1;
fast_overlap_threshold = 0.7;

% By default, only classified cells are matched
match_all = 0;

for k = 1:length(varargin)
    vararg = varargin{k};
    if ischar(vararg)
        switch lower(vararg)
            case 'noheuristic'
                use_fast_matching = 0;

            case {'matchall', 'all'}
                fprintf('%s: Matching all sources (cells and non-cells)!\n', datestr(now));
                match_all = 1;
                
                % When using both cells and non-cells, the sources are
                % expected to be denser. Let's be careful and perform an
                % exhaustive search.
                use_fast_matching = 0;
        end
    end
end

if use_fast_matching
    fprintf('%s: Using fast matching heuristic!\n', datestr(now));
else
    fprintf('%s: Fast heuristic disabled\n', datestr(now)); 
end

% Compute the matrix M of mask overlaps
%------------------------------------------------------------
num_masks = [length(masks1) length(masks2)];
M = zeros(num_masks);

available_for_match = ones(num_masks(2),1);
for i = 1:num_masks(1)
    if (mod(i,20)==0)
        fprintf('%s: Computing overlaps (%.1f%%)...\n',...
            datestr(now), 100*i/num_masks(1));
    end
    if ds1.is_cell(i) || match_all
        for j = 1:num_masks(2)
            if ds2.is_cell(j) || match_all
                if available_for_match(j)
                    M(i,j) = compute_mask_overlap(masks1{i}, masks2{j});
                    if use_fast_matching && (M(i,j) > fast_overlap_threshold)
                        available_for_match(j) = 0;
                        break;
                    end
                end
            end % j is cell
        end % j
    end % i is cell
end % i
fprintf('%s: Overlap matrix completed!\n', datestr(now));

% Find the nonzero elements of the mask overlap matrix
%------------------------------------------------------------
min_overlap_threshold = 1/3;
match_1to2 = find_nonzero_overlaps(M, min_overlap_threshold);
match_2to1 = find_nonzero_overlaps(M', min_overlap_threshold);  

end % match_masks

function match = find_nonzero_overlaps(M, threshold)

num_to_match = size(M,1);
match = cell(num_to_match,1);
for i = 1:num_to_match
    js = find(M(i,:) > threshold); % Col indices above threshold
    if isempty(js)
        match{i} = [];
    else
        overlaps = M(i,js);
        match{i} = [js' overlaps'];
        match{i} = sortrows(match{i}, 2); % Sort by overlap value
        match{i} = flipud(match{i}); % Descending
    end
end

end % find_nonzero_overlaps