function [info, masks1, masks2] = compute_affine_transform(ds1, ds2)
% IC map alignment based on user-defined control points.
%
% Inputs:
%   ds1/2: DaySummary object containing cell maps to be aligned
%
% Outputs:
%   info: Struct containing results of the affine transformation fit
%   masks1: Cell array containing thresholded masks for source 1
%   masks2: Cell array containing _transformed_ thresholded masks
%           for source 2. (NOTE: the image dimensions of masks2 match the
%           dimensions of source 1 masks!)
%

num_points_for_alignment = 4;

% To programmatically address either of the two DaySummary's
ds = cell(1,2);
ds{1} = ds1;
ds{2} = ds2;

% Display the two sets of ICs
figure;
ax1 = subplot(121);
ds1.plot_cell_boundaries('cells');
hold on;
title('Dataset 1');

ax2 = subplot(122);
ds2.plot_cell_boundaries('cells');
hold on;
title('Dataset 2');

% Allow the user to select the ICs used in matching
%------------------------------------------------------------
sel_colors = 'ybmc';
fprintf('compute_affine_transform: Please select %d cells from each dataset (in order)\n',...
    num_points_for_alignment);

selected_cells = zeros(num_points_for_alignment,2); % List of selected ICs
selected_centers = zeros(num_points_for_alignment, 2, 2); % XY position of each selected IC
num_selected = [0 0]; % Number of ICs selected from each dataset

% Loop until required number of ICs have been selected
while (~all(num_selected == num_points_for_alignment))
    click_xy = round(ginput(1)); % Get user click
    if (gca == ax1) % Axis 1 was clicked
        source_idx = 1;
    elseif (gca == ax2)
        source_idx = 2;
    end
    
    ic_idx = ds{source_idx}.get_cell_by_xy(click_xy, 'cells');
    if ~isempty(ic_idx) % Hit
        sel_idx = num_selected(source_idx) + 1;
        if (sel_idx <= num_points_for_alignment)
            boundary = ds{source_idx}.cells(ic_idx).boundary;
            fill(boundary(:,1),...
                 boundary(:,2),...
                 sel_colors(mod(sel_idx,length(sel_colors))+1));

            selected_center = mean(boundary,1);
            fprintf('  Dataset%d: Cell %d selected (at [%.1f %.1f])!\n',...
                source_idx, ic_idx,...
                selected_center(1), selected_center(2));
            
            selected_cells(sel_idx, source_idx) = ic_idx;
            num_selected(source_idx) = sel_idx;
            selected_centers(sel_idx,:,source_idx) = selected_center;
        else
            fprintf('  Dataset%d: No more cells needed!\n',...
                source_idx);
        end
    else % No hit
        fprintf('  Dataset%d: No cell detected at cursor!\n',...
            source_idx);
    end
end
fprintf('  All reference cells selected!\n');

% Transform Source2 onto Source1
%------------------------------------------------------------
tform = fitgeotrans(selected_centers(:,:,2),... % Moving points
                    selected_centers(:,:,1),... % Fixed points
                    'affine');

figure; % Pre-transform comparison
plot_boundaries_with_transform(ds1, 'b', 2, selected_cells(:,1), []);
plot_boundaries_with_transform(ds2, 'r', 1, selected_cells(:,2), []);
title('Pre-transform: Dataset1 (blue) vs. Dataset2 (red)');

figure; % Post-transform comparison
plot_boundaries_with_transform(ds1, 'b', 2, selected_cells(:,1), []);
plot_boundaries_with_transform(ds2, 'r', 1, selected_cells(:,2), tform);
title('Post-transform: Dataset1 (blue) vs. Dataset2 (red)');

% Prep output
%------------------------------------------------------------
masks1 = {ds1.cells.mask};
masks2 = {ds2.cells.mask};

mask1_ref = imref2d(size(masks1{1}));
for k = 1:ds2.num_cells
    masks2{k} = imwarp(masks2{k}, tform, 'OutputView', mask1_ref);
end

info.alignment.num_points = num_points_for_alignment;
info.alignment.selected_cells = selected_cells;
info.alignment.selected_centers = selected_centers;
info.tform = tform;

end % compute_affine_transform