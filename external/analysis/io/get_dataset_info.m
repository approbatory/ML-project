function [dataset_size, dataset_type] = get_dataset_info(hdf5_source, dataset_name)
% Determine the size and type of the HDF5 dataset

h5att = h5info(hdf5_source, dataset_name);
dataset_size = h5att.Dataspace.Size;
dataset_type = datatype_hdf5_to_matlab(h5att.Datatype.Type);