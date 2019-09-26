function spatial_weights_sparse = sparsify_spatial_weights(spatial_weights_full, max_nz_frac, verbose)

if ~exist('max_nz_frac', 'var') || isempty(max_nz_frac)
    max_nz_frac = nnz(spatial_weights_full)/numel(spatial_weights_full);
end

if ~exist('verbose', 'var')
    verbose = true;
end

spatial_weights_sparse = ndSparse.spalloc(size(spatial_weights_full),...
    numel(spatial_weights_full)*max_nz_frac);

for i = 1:size(spatial_weights_full, 3)
    spatial_weights_sparse(:,:,i) = ndSparse(spatial_weights_full(:,:,i));
    if verbose
        fprintf('Done %d / %d\n', i, size(spatial_weights_full, 3));
    end
end

end
