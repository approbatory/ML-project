function setup

addpath('utils', 'plotting', 'decoding', 'data_loading', 'supfigs_code', 'testing_pf_width_discrepancy');
addpath(genpath('external'));
addpath .

%if ispc
%    DecodeTensor.linear_track_path('../../../Box/share');
%end