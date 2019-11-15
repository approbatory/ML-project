function ax = ClearFigure;
% ax = ClearFigure;
%
% Clears the current figure, removes the menu and tool bars, colors
% the background white, and adds a inserts a very large set of axes.
clf;
set(gcf,'MenuBar','none','ToolBar','none','Color',[1 1 1]);
%Q = ComputeSubplotPositions(1,1,{},0.05,0.01,0,0.05,0.05, 0);
%ax = subplot('Position',Q);