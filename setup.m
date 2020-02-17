function setup(p)

addpath('utils', 'plotting', 'decoding', 'data',_'loading'
addpath(genpath('external'));
addpath .

if ispc
    DecodeTensor.linear_track_path('../../../Box/share');
end