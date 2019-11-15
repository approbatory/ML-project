function TightenAxesToFigure
% TightenAxesToFigure
%
% Tightens the current axes to the bounds of the figure by
% repositioning the axes so that there are no margins on the sides.
% 
% This is done by using the undocumented axis property
% 'LooseInset'. To avoid having to call this function everytime the
% figure is resized, ResizeFcn of the figure is set to update
% LooseInset in the same way.

set(gca,'LooseInset',get(gca,'TightInset'))
set(gcf,'ResizeFcn',@(varargin) set(gca,'LooseInset',get(gca,'TightInset')));

