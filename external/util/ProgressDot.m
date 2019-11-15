function ProgressDot(n, numTotal, numPerLine)
% ProgressDot(n, numTotal, numPerLine)
%
% Prints a dot. If n == 0 (mod numPerLine) | n == numTotal, prints
% 'n' as well and then a newline.
fprintf('|');
if (mod(n, numPerLine) == 0 | n == numTotal)
  fprintf('%d/%d\n', n, numTotal);
end