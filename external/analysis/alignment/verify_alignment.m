function verify_alignment(ds1, ds2, affine_info, m_1to2, m_2to1)

ds1_matched = find(~cellfun(@isempty, m_1to2));
ds2_matched = find(~cellfun(@isempty, m_2to1));

subplot(121);
plot_boundaries_with_transform(ds2, 'r', 1, [], affine_info.tform);
plot_boundaries_with_transform(ds1, 'b', 1, ds1_matched, []);
axis equal tight;
title('Matched cells from Dataset1 (blue fill)');

subplot(122);
plot_boundaries_with_transform(ds1, 'b', 1, [], []);
plot_boundaries_with_transform(ds2, 'r', 1, ds2_matched, affine_info.tform);
axis equal tight;
title('Matched cells from Dataset2 (red fill)');