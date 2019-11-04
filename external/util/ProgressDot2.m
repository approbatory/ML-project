function ret = ProgressDot2(n, numTotal, numPerDot, dotsPerLine, varargin)
% ret = ProgressDot(n, numTotal, numPerDot, dotsPerLine, varargin)
%
% Prints a dot. Every 'numPerDot' iterations it prints a dot, and after
% dotsPerLine dots it prints the current count and a new line.

numPerLine = numPerDot*dotsPerLine;

if (mod(n, numPerDot) == 0)
  fprintf('|');
end

ret =  0;
if (mod(n, numPerLine) == 0 | n == numTotal)
  suffix = [];
  if (~isempty(varargin) & ~isempty(varargin{1}))
    suffix = varargin{1}(varargin{2:end});
  end
  fprintf('%d/%d %s\n', n, numTotal, suffix);
  ret = 1;
end