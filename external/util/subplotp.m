function h = subplotp(Q,n)
% h = subplotp(Q,n)
%
% Creates a subplot using the position information in the n'th row
% of Q, and returns a handle to it. Useful when combined with
% ComputeSubplotPositions.
%
% See also: ComputeSubplotPositions.

if (isempty(Q))
  error('Position matrix is empty.');
end

if (size(Q,2) ~= 4)
  error('Position matrix should have four columns.');
end

if (n<1 | n> size(Q,1))
  error('Row index into position matrix should be in the range [1, %d].', size(Q,1));
end

h = subplot('Position',Q(n,:));