function next_idx = find_next_cell_to_process(current_idx, to_process)
% Helper function to identify the next cell to process, used by
% `classify_cells` and others.
%
% `to_process` is a Boolean vector of length equal to the number of cells,
% where 0 indicates that the cell has already been processed and 1
% indicates that the cell still needs to be processed.
% 
% The function returns the first cell index larger than 'current_idx' that
% has yet to be processed, and loops at the end. If all cells have been
% processed, then returns an empty matrix.

if any(to_process)
    num_cells = length(to_process);

    to_process = circshift(to_process, -current_idx);
    next_offset = find(to_process, 1);
    next_idx = mod(current_idx + next_offset - 1, num_cells) + 1;
else
    next_idx = [];
end
