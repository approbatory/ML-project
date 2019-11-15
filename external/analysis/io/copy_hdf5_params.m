function copy_hdf5_params(hdf5_from, hdf5_to)
% Transfer all datasets under the HDF5 directory '/Params' from one HDF5
% file to another.
%
% Inputs:
%   - hdf5_from: Name of the source HDF5 file
%   - hdf5_to:   Name of the target HDF5 file
%

% Enumerate all 'Params' datasets in the source file
params_root = '/Params';
try
    params = h5info(hdf5_from, params_root);
catch
    fprintf('Warning: Source HDF5 file lacks "Params" directory\n');
    return;
end
params = params.Datasets;
num_params = length(params);

for i = 1:num_params
    param = params(i);
    param_path = sprintf('%s/%s', params_root, param.Name);
    param_type = datatype_hdf5_to_matlab(param.Datatype.Type);
    
    % Read the parameter value from source
    param_val = h5read(hdf5_from, param_path);
    
    % Copy the parameter value to target
    try
        h5create(hdf5_to,...
                param_path,...
                param.Dataspace.Size,...
                'Datatype', param_type);
        h5write(hdf5_to, param_path, cast(param_val, param_type));
    catch e
        fprintf('  Warning: %s\n', e.message);
    end
end