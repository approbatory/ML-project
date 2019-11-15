function lsm
% lsm
%
% Lists the MATLAB files in the current directory. First the m files,
% follows by the mat files, and finally the mex files.

disp('MATLAB source files:');
try
  ls -ltr *.m
catch
  disp('  None found.');
end

disp('MATLAB data files:');
try
  ls -ltr *.mat
catch
  disp(' None found.');
  disp(' ');
end

disp('MATLAB mex files:');
try
  ls -ltr *.mex*
catch
  disp('  None found.');
  disp(' ');
end
