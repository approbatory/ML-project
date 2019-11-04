function [hl, hp] = PlotMeanWithFilledErrorBand(x, y, pe, ne, lineColor, lineWidth, bandColor,bandAlpha)
% PlotMeanWithFilledErrorBand(x, y, pe, ne, lineColor, lineWidth, bandColor, bandAlpha)
%
% Plots (x,y) with the specified line color and line width, on top of
% a filled band whose extent at x(i) is [y(i)+pe(i), y-ne(i)]. The
% band has the specified color and alpha.
%
% Example: Plot mean and std of 20 trials of 100 dimensional random variable drawn from N(0,1^100).
%
% r = randn(100,20);
% m = mean(r,2);
% s = std(r,[],2);
%
% PlotMeanWithFilledErrorBand(1:100,m,s,s,'r',2,'r',0.25);
x = x(:)';
y = y(:)';
pe = pe(:)';
ne = ne(:)';

ub = y+pe;
lb = y-ne;

u = [x(1:end-1); x(2:end); x(2:end); x(1:end-1)];
v = [lb(1:end-1); lb(2:end); ub(2:end); ub(1:end-1)];
hp = patch(u,v,bandColor,'FaceAlpha', bandAlpha, 'EdgeColor','none','LineStyle','none');
hold on;
hl = plot(x,y,'LineWidth', lineWidth, 'Color', lineColor);


