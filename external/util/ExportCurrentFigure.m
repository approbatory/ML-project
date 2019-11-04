function ExportCurrentFigure(figureId, figureName, varargin)
% ExportCurrentFigure(figureId, figureName, varargin)
%
% Exports the current figure to pdf file with name
% figFIGUREID_FIGURENAME. An optional argument specifies the desired
% output format (e.g. pdf, eps, etc.).
%
% Example:
%
% >> clf
% >> plot(randn(1000));
% >> ExportCurrentFigure(1,'randnExample', 'pdf');
% Exported "fig1_randnExample.pdf".
% 
% See also: export_fig.

if (isempty(varargin))
  outputFormat = 'pdf';
else
  outputFormat = varargin{1};
end
if (isnumeric(figureId))
  fullName = sprintf('fig%d_%s', figureId, figureName);
else
  fullName = sprintf('fig%s_%s', figureId, figureName);
end
fullName = strrep(fullName,' ','_');
eval(sprintf('! rm -f %s.*', fullName));
eval(sprintf('export_fig ''%s'' -%s', fullName, outputFormat));
fprintf('Exported "%s.%s".\n', fullName, outputFormat);