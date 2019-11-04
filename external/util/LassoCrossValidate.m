function [indBestLambda, mu, sem] = LassoCrossValidate(X,y,lambda,k,doPlot, varargin)
% [indBestLambda, mu, sem] = LassoCrossValidate(X,y,lambda,k,doPlot, varargin)
%
% Cross-validation function for LassoRegress. Performs k-fold cross
% validation on the regressors X and dependent variable for the
% penalties in lambda, and returns the mean and sd of sse for each
% lambda value. If DOPLOT = 1, a plot of the values is also
% made. Additional arguments are passed onto LassoRegress;

p = inputParser;
p.KeepUnmatched = true; % Additional varargs will be passed onto LASSOF, so don't complain if you get some you don't expect.
p.addOptional('LassoFunction',@LassoRegress);
p.parse(varargin{:});

LASSOF = p.Results.LassoFunction;

numSamples = size(X,1);
cvp = cvpartition(numSamples,'k',k);

mu = zeros(1,k);
sd = zeros(1,k);

sse = zeros(k,numel(lambda));
parfor i = 1:k
  itrn = cvp.training(i);
  itst = cvp.test(i);
  [B,b0] = LASSOF(X(itrn,:),y(itrn),'lambda',lambda,varargin{:});
  ytst = y(itst);
  Yfit = bsxfun(@plus, X(itst,:)*B,b0);
  Yerr = bsxfun(@minus,y(itst),Yfit);  
  sse(i,:) = sum(Yerr.^2,1);
end

mu = mean(sse,1);
sd = std(sse,[],1);
sem = sd/sqrt(k);

indBestLambda = max(argmin(mu)); % In case there are multiple mins, take the one with the highest constraint.

if (doPlot)
  ttl = 'Lasso Cross Validation Results';
  % Check to see if the figure exists
  hfig = findobj('name',ttl);
  if (isempty(hfig))
    hfig = sfigure;
  end
  set(0,'CurrentFigure',hfig);
  set(hfig,'Name',ttl);
  clf;
  u = lambda(:)';
  v = mu(:)';
  s = sem(:)';
  loglog(u,v,'bo','MarkerSize',5,'MarkerFaceColor','b'); hold on;
  line([u;u],[v-3*s;v+3*s],'Color','b');
  loglog(u,v,'b');
  loglog(u(indBestLambda),v(indBestLambda),'ro','MarkerFaceColor','r','MarkerSize',5);
end


