addpath utils plotting decoding data_loading
addpath(genpath('external'));

if ispc
    DecodeTensor.linear_track_path('../../../Box/share');
end