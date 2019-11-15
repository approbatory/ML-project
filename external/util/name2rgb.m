function rgb = name2rgb(name,varargin)
% rgb = name2rgb(name,[colorCodesFile])
%
% Attempts to convert the specified (case-insensitive) name into an
% rgb value using a lookup table of 500 named colors, taken from
% http://cloford.com/resources/colours/500col.htm
%
% See also: COLORSWATCH.

if (isempty(varargin))
  colorCodesFile = fullfile(fileparts(mfilename('fullpath')), 'colorCodes.mat');
else
  colorCodesFile = varargin{1};
end

load(colorCodesFile);
  
rgb = [];
name = lower(name);
for i = 1:numel(colorCodes)
  if (isequal(name, colorCodes(i).name))
    rgb = colorCodes(i).rgb/255;
  end
end
if (isempty(rgb))
  error('Unrecognized color name "%s".', name);
end
    