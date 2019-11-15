function Y = when(varargin)
% Y = when(p1, v1, p2, v2,...,velse)
%
% If p1, then Y = v1, elseif p2, Y = v2,..., else Y = velse.
%
% Note that MATLAB evaluate the arguments get evaluate even if they
% don't end up being used.

if (nargin<3)
  error('Expected atleast three arguments.');
elseif(mod(nargin,2)~=1)
  error('Expected an odd number of arguments >= 3.');
end

assignedY = 0;
for i = 1:2:numel(varargin)-1
  if (varargin{i})
    Y = varargin{i+1};
    assignedY = 1;
    break;
  end
end

if (~assignedY)
  Y = varargin{end};
end

