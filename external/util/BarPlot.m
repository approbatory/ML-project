function [hb, hl] = BarPlot(c, lb, ub, hw, colors,edges)
% [hb, hl] = BarPlot(c, lb, ub, hw, colors,edges)
%
% Produces a bar plot in the current axes by plotting, for the i'th
% element of each of the input vectors, a bar of height c(i) centered
% at i, with width 2*hw, and of color colors(i,:). A black line is
% also drawn from lb(i) to ub(i). If edges is non-zero, edges are
% drawn black, otherwise are the color of the bar.
% 
% HB is a vector of handles to the bar plots themselves, and HL is a
% vector of handles to the lines.
hb = [];
hl = [];
for i=1:numel(c)
  hb(i) = patch([i-hw;i-hw;i+hw;i+hw],[0 c(i) c(i) 0],colors(i,:));
  if (~edges)
    set(hb(i),'EdgeColor',colors(i,:));
  end
  hold on;
end

% Plot the lines in their own loop to avoid messing up subsequent
% legends.
hl = [];
for i = 1:numel(c)
  hl(i) = line([i,i],[lb(i) ub(i)],'Color','k');
end
