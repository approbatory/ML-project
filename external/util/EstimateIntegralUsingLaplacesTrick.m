function [Q, xpk, fpk, H, exitFlag] = EstimateIntegralUsingLaplacesTrick(f, x0, varargin)
% [Q, xpk, fpk, H, exitFlag] = EstimateIntegralUsingLaplacesTrick(f, x0, varargin)
%
% Estimates the integral of the function F by using Laplace's trick to
% approximate it by a multidimensional normal distribution centered at
% its peak, with covariance matrix equal to the inverse of the hessian
% of the log of F at the peak. The estimate is returned in Q. The peak
% is found numerically with starting point X0. XPK contains the
% coordinates of the peak of F, FPK is its value there. H is the
% hessian at the peak, with the peak normalized to 1.  If an optional
% input argument is set to 'logf', then f is assumed to be the
% logarithm of the function, not the function itself.

if (numel(varargin)>0 && isequal(lower(varargin{1}),'logf'))
  lf = f;
  f = @(x) exp(lf(x));
else
  lf = @(x) log(f(x));
end
nd = numel(x0); % number of dimensions

warning off;
opts = optimset('MaxFunEvals',100000,'MaxIter',100000);
[xpk,lfmin,exitFlag] = fminsearch(@(x) -lf(x), x0, opts); % Find the peak.
fpk = f(xpk);  % The value at the peak.
H   = hessian(lf, xpk); % Note that we're taking the hessian of the log of the function
H(isnan(H)) = 0;
if (exitFlag>0)
  Q = fpk*(sqrt(2*pi)^nd)/sqrt(abs(det(H)));
else
  Q = 0;
end





