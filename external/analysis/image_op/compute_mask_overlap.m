function overlap = compute_mask_overlap(mask1, mask2)
% Compute the Jacard similarity of two (logical) masks. The two masks
%   should have the same dimensions.

overlap = nnz(mask1 & mask2) / nnz(mask1 | mask2);