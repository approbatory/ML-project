function PlotCustom3DGrid(ax, lineStyle,plotBoxEdge)
% PlotCustom3DGrid(ax, lineStyle,[plotBoxEdge])
%
% Plots a grid in the specified axes with the specified
% linestyle. Setting plotBoxEdge to non-zero will also plot edges
% around the frame, in black.
%
% Example: Make the current grid red and dotted
%
% PlotCustum3DGrid(gca,'LineStyle',':','Color','r');

xtk = get(ax,'xtick');
ytk = get(ax,'ytick');
ztk = get(ax,'ztick');
% First coordinate is moving, second is anchoring
hxy = line(xlim'*ones(size(ytk)),[ytk;ytk],min(zlim)*ones(size([ytk;ytk;])));
hyx = line([xtk;xtk],ylim'*ones(size(xtk)),min(zlim)*ones(size([xtk;xtk;])));
hzy = line(min(xlim)*ones(size([ytk;ytk])),[ytk;ytk],zlim'*ones(size(ytk)));
hyz = line(min(xlim)*ones(size([ztk;ztk])),ylim'*ones(size(ztk)),[ztk;ztk]);
hxz = line(xlim'*ones(size(ztk)),min(ylim)*ones(size([ztk;ztk])),[ztk;ztk]);
hzx = line([xtk;xtk],min(ylim)*ones(size([xtk;xtk])),zlim'*ones(size(xtk)));

h = {hxy,hyx,hzy,hyz,hxz,hzx};
cellfun(@(h) set(h,lineStyle{:}), h);

if (nargin>2)
  if (plotBoxEdge)
    line(xlim,min(ylim)*[1 1], max(zlim)*[1 1],'Color','k');
    line(min(xlim)*[1 1],ylim, max(zlim)*[1 1],'Color','k');
    line(min(xlim)*[1 1],max(ylim)*[1 1], zlim,'Color','k');
  end
end