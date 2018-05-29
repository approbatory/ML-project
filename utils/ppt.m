function ppt(XY, varargin)
plot3(1:size(XY,1),XY(:,1), XY(:,2), varargin{:});