function H = serrorbar(varargin)
assert(nargin ~= 0, 'missing input');
color = [];
if ischar(varargin{end})
    color = varargin{end};
    varargin{end} = [];
    nargs = nargin - 1;
else
    nargs = nargin;
end


if nargs == 2
    x = 1:length(varargin{1});
    y = varargin{1};
    e = varargin{2};
    if isempty(color)
        H = shadedErrorBar(x,y,e);
    else
        H = shadedErrorBar(x,y,e,'lineProps',color);
    end
elseif nargs == 3
    x = varargin{1};
    y = varargin{2};
    e = varargin{3};
    if isempty(color)
        H = shadedErrorBar(x,y,e);
    else
        H = shadedErrorBar(x,y,e,'lineProps',color);
    end
else
    x = varargin{1};
    y = varargin{2};
    e = varargin{3};
    extra = varargin(4:end);
    if isempty(color)
        H = shadedErrorBar(x,y,e,extra{:});
    else
        H = shadedErrorBar(x,y,e,'lineProps',color,extra{:});
    end
end

