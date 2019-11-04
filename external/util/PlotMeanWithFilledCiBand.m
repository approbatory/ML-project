function [hl, hb] = PlotMeanWithFilledCiBand(x, y, yu, yl, lineColor, lineWidth, bandColor,bandAlpha)
% [hl, hb] = PlotMeanWithFilledCiBand(x, y, yu, yl, lineColor, lineWidth, bandColor, bandAlpha)
%
% Plots (x,y) with the specified line color and line width, on top of
% a filled band whose extent at x(i) is [yl(i), yu(i)]. The band has
% the specified color and transparency alpha. This is useful for
% plotting confidence intervals as a band.
%
% HL is a handle to the line, and HB is a handle to the patch forming
% the band.
%
% Example: Plot mean and std of 20 trials of 100 dimensional random
% variable drawn from N(0,1^100).
%
% r = randn(100,20);
% m = mean(r,2);
% pu = prctile(r,95,2);
% pl = prctile(r,5,2);
%
% figure(1);
% clf;
% PlotMeanWithFilledCiBand(1:100,m,pu,pl,'r',2,'r',0.25);
%
% See also: PLOTMEANWITHFILLEDERRORBAND.

x = x(:)';
y = y(:)';
yu = yu(:)';
yl = yl(:)';

u = [x(1:end-1); x(2:end); x(2:end); x(1:end-1)];
v = [yl(1:end-1); yl(2:end); yu(2:end); yu(1:end-1)];
hb = patch(u,v,bandColor,'FaceAlpha', bandAlpha, 'EdgeColor',bandColor,'EdgeAlpha',bandAlpha,'LineStyle','none');
hold on;
hl = plot(x,y,'LineWidth', lineWidth, 'Color', lineColor);

