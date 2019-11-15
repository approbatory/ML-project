function [match_1to2, match_2to1] = bijective_filter(match_1to2, match_2to1)
% Makes sure that the the matching is bijective, i.e.
%   1. Each cell in DatasetI matches to at most one cell in DatasetJ.
%   2. Matches are kept only when they are bidirectional.
%
% Assumes that each element of match_itoj are descending in their overlap
% measure (see 'match_masks.m').
%

% Step 1: Each cell can keep only its highest matched partner
%------------------------------------------------------------
num_cells1 = length(match_1to2);
num_cells2 = length(match_2to1);

for i = 1:num_cells1
    match_itoj = match_1to2{i};
    num_matches = size(match_itoj,1);
    if (num_matches > 1)
        match_1to2{i} = match_itoj(1,:);
    end
end

for j = 1:num_cells2
    match_jtoi = match_2to1{j};
    num_matches = size(match_jtoi,1);
    if (num_matches > 1)
        match_2to1{j} = match_jtoi(1,:);
    end
end

% Step 2: We now have a bipartite directed graph, with at most one edge
% from each node. Next, we filter for pairings where the connections are
% bidirectional, i.e. keep pairings between Cell X in Dataset1 and Cell Y
% in Dataset2 if and only if Cell Y is the highest overlap match for Cell X
% and Cell X is the highest overlap match for Cell Y.
%------------------------------------------------------------
for i = 1:num_cells1
    match_itoj = match_1to2{i};
    if ~isempty(match_itoj)
        % Check if the directed edge is reciprocated
        j = match_itoj(1);
        match_jtoi = match_2to1{j};
        if ~isempty(match_jtoi)
            if (match_jtoi(1) ~= i) % Unreciprocated
                match_1to2{i} = [];
            end
        else
            match_1to2{i} = [];
        end
    end
end

for j = 1:num_cells2
    match_jtoi = match_2to1{j};
    if ~isempty(match_jtoi)
        i = match_jtoi(1);
        match_itoj = match_1to2{i};
        if ~isempty(match_itoj)
            if (match_itoj(1) ~= j)
                match_2to1{j} = [];
            end
        else
            match_2to1{j} = [];
        end
    end
end

end % bijective_filter