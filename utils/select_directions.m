function [fw, bw, cm_per_pix] = select_directions(XY, varargin)
p = inputParser;
p.addRequired('XY', @(x) isnumeric(x) && size(x,2) < 3);
p.addParameter('speed_threshold', 4, @isnumeric);
p.addParameter('total_length', 118, @isnumeric);
p.addParameter('frame_sampling_freq', 20, @isnumeric);
p.parse(XY, varargin{:});


track_coord = p.Results.XY(:,1);
cm_per_pix = p.Results.total_length ./ range(track_coord);
velocity = [0; diff(track_coord)] .*...
    cm_per_pix .*...
    p.Results.frame_sampling_freq;
fw = velocity > p.Results.speed_threshold;
bw = velocity < -p.Results.speed_threshold;
end