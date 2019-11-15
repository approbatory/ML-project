function h = FindFigureCreate(name)
% function h = FindFigureCreate(name)
%
% Returns the handle to the figure with the specified name, or creates
% one if it can't be found. The search is done by tags not names, so
% the name of the figure can be changed afterwards.

h = findobj('tag',name);
if (isempty(h))
  h = figure;
  set(h,'name',name,'tag',name);
elseif (numel(h)>1)
  warning('Multiple figures found with tag "%s".', name);
  h = h(1);
end
