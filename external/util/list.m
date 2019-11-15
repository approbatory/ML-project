function list(varargin)
% list(...)
%
% If one argument is provided, evaluates
%
% find ~ -name "arg1"
%
% If two arguments are provided, evaluates
%
% find arg1 -name "arg2"

switch(numel(varargin))
 case 1
  eval(sprintf('! find ~ -name "%s" | grep -v "todelete" | grep -v "research_backup"', varargin{1}));
 case 2
  eval(sprintf('! find %s -name "%s"| grep -v "todelete" | grep -v "research_backup"', varargin{1}, varargin{2}));
 otherwise
  error('Expected 1 or 2 arguments.');
end

