function [match_1to2, match_2to1, info] = run_alignment(ds1, ds2, varargin)
% Align two sets of cell filters.
%
% Inputs:
%   ds1/2: DaySummary object containing cell maps to be aligned
%
% Optional inputs:
%   'notrans': No transformation applied between ds1 and ds2
%   'keepall': A cell from ds1 may match to more than one cell in ds2, and
%              visa versa.
%
%   Additional optional arguments are passed into 'match_masks.m', e.g.
%       'fast' and 'matchall'
%
% Outputs:
%   match_XtoY: Cell that contains mapping information from source X to
%       source Y.
%
%       For example, match_1to2{k} is a [Nx2] matrix where N is the number
%       of ICs from source 2 that match cell k of source 1. The first column
%       of the matrix is the index of the matching cell in source 2; the
%       second column is the overlap score between the cell pairs.
%
%   info: Struct with additional information (e.g. selected points, etc.)
%       regarding the alignment run.
%
% Example usage:
%   [m1to2, m2to1] = run_alignment(m1d12, m1d13, 'fast')
%

% Default alignment options
%------------------------------------------------------------
wait_for_prompt = 1;
use_transform = 1;
bijective_matching = 1;

for k = 1:length(varargin)
    if ischar(varargin{k})
        switch varargin{k}
            case 'notrans'
                % No affine transform needed (e.g. for matching multiple
                % extraction runs on the same movie). Requires the image
                % dimensions to be identical!
                assert(all(size(ds1.cells(1).mask)==size(ds2.cells(1).mask)),...
                    'notrans option requires cell image dimensions to be identical!');
                use_transform = 0;
            case 'keepall' % Keep all matches
                bijective_matching = 0;
            case 'noprompt' % For programmatic use
                wait_for_prompt = 0;
        end
    end
end

% If the DaySummary has no explicitly classified cells, assume the user 
% wants to match all cells
ds1_set_labels = all(cellfun(@isempty, {ds1.cells.label}));
ds2_set_labels = all(cellfun(@isempty, {ds2.cells.label}));
if ds1_set_labels
    ds1.set_labels;
end
if ds2_set_labels
    ds2.set_labels;
end

% Control point-based registration of two sets of ICs
%------------------------------------------------------------
if use_transform
    fprintf('run_alignment: Beginning alignment...\n');
    [info, masks1, masks2] = compute_affine_transform(ds1, ds2);
else
    masks1 = {ds1.cells.mask};
    masks2 = {ds2.cells.mask};
    info = [];
    
    figure;
    plot_boundaries_with_transform(ds1, 'b', 2, [], []);
    plot_boundaries_with_transform(ds2, 'r', 1, [], []);
    title('Dataset1 (blue) vs. Dataset2 (red)');
end

if wait_for_prompt
    input('run_alignment: Press enter to continue with mask matching >> ');
end

[match_1to2, match_2to1] = match_masks(masks1, masks2, ds1, ds2, varargin{:});

% Optional bijective filtering
%------------------------------------------------------------
if bijective_matching
    fprintf('run_alignment: Applying bijective filter...\n');
    [match_1to2, match_2to1] = bijective_filter(match_1to2, match_2to1);
end

info.num_matches = sum(~cellfun(@isempty, match_1to2));
fprintf('run_alignment: Found %d matches!\n', info.num_matches);

% "Transfer" unmatched filters across datasets. TODOs:
% - Clean up code
% - Check if the transferred filter is out of bounds
%------------------------------------------------------------
if use_transform % Transfer is irrelevant for in-place matches
    % ds2 --> ds1
    unmatched_from_ds2 = find(cellfun(@isempty, match_2to1)' & ds2.is_cell);
    num_unmatched_from_ds2 = length(unmatched_from_ds2);
    
    ds1_imsize = size(ds1.cells(1).im);
    ds1_ref = imref2d(ds1_imsize);
    filters_2to1 = zeros(ds1_imsize(1), ds1_imsize(2), num_unmatched_from_ds2);
    coms_2to1 = zeros(num_unmatched_from_ds2, 2);
    for k = 1:num_unmatched_from_ds2
        cell_idx2 = unmatched_from_ds2(k);
        filters_2to1(:,:,k) = imwarp(...
            ds2.cells(cell_idx2).im, info.tform, 'OutputView', ds1_ref);
        coms_2to1(k,:) = transformPointsForward(info.tform, ds2.cells(cell_idx2).com');
    end
    info.filters_2to1.ds2_inds = unmatched_from_ds2;
    info.filters_2to1.im = filters_2to1;
    info.filters_2to1.com = coms_2to1;
    
    % ds1 --> ds2
    unmatched_from_ds1 = find(cellfun(@isempty, match_1to2)' & ds1.is_cell);
    num_unmatched_from_ds1 = length(unmatched_from_ds1);
    
    ds2_imsize = size(ds2.cells(1).im);
    ds2_ref = imref2d(ds2_imsize);
    invtform = invert(info.tform); % tform is computed for ds2 --> ds1
    filters_1to2 = zeros(ds2_imsize(1), ds2_imsize(2), num_unmatched_from_ds1);
    coms_1to2 = zeros(num_unmatched_from_ds1, 2);
    for k = 1:num_unmatched_from_ds1
        cell_idx1 = unmatched_from_ds1(k);
        filters_1to2(:,:,k) = imwarp(...
            ds1.cells(cell_idx1).im, invtform, 'OutputView', ds2_ref);
        coms_1to2(k,:) = transformPointsForward(invtform, ds1.cells(cell_idx1).com');
    end
    info.filters_1to2.ds1_inds = unmatched_from_ds1;
    info.filters_1to2.im = filters_1to2;
    info.filters_1to2.com = coms_1to2;
end

% If 'run_alignment' internally set the DaySummary labels for matching,
% then undo before exiting.
if ds1_set_labels
    ds1.reset_labels;
end
if ds2_set_labels
    ds2.reset_labels;
end