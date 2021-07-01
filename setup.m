function setup

addpath('utils', 'plotting', 'decoding', 'data_loading', 'supfigs_code');
addpath(genpath('external'));
addpath .

%if ispc
%    DecodeTensor.linear_track_path('../../../Box/share');
%end