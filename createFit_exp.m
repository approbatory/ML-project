function [fitresult, gof] = createFit_exp(n, m)
[xData, yData] = prepareCurveData(n, m);
ft = fittype('I_0*N*(1 - 2.^(-x/N))', 'independent', 'x', 'dependent', 'y');
opts = fitoptions('Method', 'NonlinearLeastSquares');
opts.Display = 'Off';
opts.StartPoint = [1 100];
[fitresult, gof] = fit(xData, yData, ft, opts);