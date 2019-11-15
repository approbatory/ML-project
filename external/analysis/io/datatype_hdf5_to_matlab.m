function datatype_matlab = datatype_hdf5_to_matlab(datatype_hdf5)
% Convert the HDF5 datatype names to Matlab equivalents

datatype_matlab = '';
switch (datatype_hdf5)
    case 'H5T_STD_I16LE'
        datatype_matlab = 'int16';
    case 'H5T_STD_U16LE'
        datatype_matlab = 'uint16';
    case 'H5T_IEEE_F64LE'
        datatype_matlab = 'double';
    case 'H5T_IEEE_F32LE'
        datatype_matlab = 'single';
end