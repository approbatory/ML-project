function [B, b0, lambda, stats] = LassoRegress2NonNegative(X,y,varargin)
% [B, b0, lambda, stats] = LassoRegress2NonNegative(X,y,varargin)
%
% A faster version of LassoRegress, that uses Covariance Updates, and
% requires the weights to be positive.
% 
% Computes the dual function of the primal Lasso regression problem:
%
% minimize    (y - b0 - Xb)'(y - b0 - Xb)
% subject to  |b|_1 <= t
%
% The Lagrangian of this function is
%
% L(b0,b,s) = (y - b0 - Xb)'(y - b0 - Xb) + s(|b|_1 - t)
%    s >= 0
%
% Hence the Lagrangian dual function is
%
% g(s) = inf_(b0,b)  (y - b0 - Xb)'(y - b0 - Xb) + s(|b|_1 - t)
%
% Here we try to estimate the values of b0, b that minimize g(s), for
% a given value of s. Since g(s) is a lower bound for the primal
% problem, the idea is then to find the value of s that maximizes
% g(s), although it seems that what is typically done instead is to
% use cross validation to find the s whos (b0,b) values produce the
% lowest cv-error.
%
% The algorithm works as follows (based on [2], and on MATLAB's 'lasso' function). 
%
% 1. The predictors and the dependent variable are first centered.
%
% 2. We initially start with no active variables.
%
% 3. In the main loop, cyclic coordinate descent is applied to the
%    active set.  This is continued until convergence of the weights.
%
% 4. Once weight convergence is achieved, we check for potentially
%    active weights, add them to the list, and run coordinate descent
%    again.
%
% 5. If the active set doesn't change, or doesn't change by enough,
%    the algorithm terminates, other wise it continues.
%
% The weights are returned in the vector b, the constant offsets in
% b0, and additional statistics are returned in the structure STATS.
%
% References:
%
% [1] Boyd, "Convex Optimization"
%
% [2] Friedman et. al "Regularization Paths for GLMs via Coordinate
% Descent", Journal of Statistical Software 2007
%
% See also: LassoRegress, LassoRegress2.

Args = inputParser;
Args.KeepUnmatched = true;
Args.addOptional('lambda',[]);
Args.addOptional('predictorWeights',[]);
Args.addOptional('initial',[]);
Args.addOptional('numLambda',100);
Args.addOptional('lambdaRatio',1e-4);
Args.addOptional('relTol',1e-4);
Args.addOptional('verbose','no');
Args.parse(varargin{:});
Args = Args.Results;

relTol = Args.relTol;
verbose = isequal(lower(Args.verbose),'yes');

numRegressors = size(X,2);

predictorWeights = abs(Args.predictorWeights);
if (isempty(predictorWeights))
  predictorWeights = ones(numRegressors,1);
end

ym = mean(y);
Xm = mean(X);

X  = bsxfun(@minus,X,Xm);
Xss = sum(X.^2,1);
y  = y - ym;

lambda = Args.lambda;
if (isempty(lambda))
  numLambda = Args.numLambda;
  if (numLambda<=0)
    error('No lambda value provided an numLambda is not positive.');
  end
  % Compute the maximum lambda value.  At this value, only one
  % regressor will be allowed in the model.  The magnitude of a
  % (non-penalized) update to this regressor will be 2*r'*X, where r
  % is the residuals without this regressor. But since it's the only
  % one on the model, r = y (since we mean subtracted).
  lambdaMax = max(abs(2*y'*X(:,predictorWeights~=0))); 
  lambdaMin = Args.lambdaRatio*lambdaMax;
  lambda = logspace(log10(lambdaMin),log10(lambdaMax),numLambda);
end

lambda = sort(lambda,'descend');
numLambda = numel(lambda);

B = zeros(numRegressors, numLambda);
b = zeros(numRegressors,1); % This is where the coefficients will be stored for a given lambda.
active = false(1,numRegressors);

initial = Args.initial;
useInitial = false;
if (isempty(initial))
  useInitial = false;
else
  if (size(initial,2) == numLambda)
    useInitial = true;
    initialVals = initial;
  elseif (size(initial,2)==1) % Single vector supplied, use it for all
    useInitial = true;
    initialVals = initial*ones(1,numLambda);    
  else
    error('Initial value supplied but size ~= 1 or the number of lambda values supplied.');
  end
end

if (useInitial && ~isequal(size(initialVals),[numRegressors, numLambda]))
  error('Size mismatch for initial values.');
end

iters = ones(numLambda,1);
df = 0*iters;
sse = 0*iters;

Xy = y'*X;
XX = X'*X;

for i = 1:length(lambda)
  if (verbose)
    fprintf('Running for lambda = %1.3f\n', lambda(i));
  end

  if (useInitial)
    b = initialVals(:,i);
  else
    if (i>1) % Use a warm start
      b = bnext;
      active = activeNext;
    end
  end

  while(1)
    [b1,active1] = cyclicCoordinateDescent(X,y,b,Xss,active,lambda(i),predictorWeights,XX,Xy);
    if (verbose)
      fprintf(' Iteration %d: SSE = %1.3f\n', iters(i), sum((y - X*b).^2));
    end
    % Check for convergence
    maxRelChange = max(abs(b1-b)./max(max(abs(b),abs(b1)),1));
    if (verbose)
      fprintf(' Largest relative change: %1.5f\n', maxRelChange);
    end
    converged = maxRelChange<=relTol;
    if (converged)
      % Try one more pass
      potentiallyActive = thresholdActive(X,y,b1,active1,lambda(i),predictorWeights) | active1;
      [b2,active2] = cyclicCoordinateDescent(X,y,b1,Xss,potentiallyActive,lambda(i),predictorWeights,XX,Xy);      
      if (isequal(active2, active1)) % No change to the model, so we're done.
        bnext = b1;
        activeNext = active1;
        break; 
      elseif max(abs(b2-b1)./max(max(abs(b2),abs(b1)),1)) < relTol
        % Model didn't change by much, take the smaller model and exit.
        bnext = when(sum(active1)<=sum(active2), b1, b2);
        activeNext = when(sum(active1)<=sum(active2), active1, active2);
        break;
      else
        % Model changed significantly
        bnext = b2;
        activeNext = active2;
      end
    else
        bnext = b1;
        activeNext = active1;
    end
    b = bnext;
    active = activeNext;
    iters(i) = iters(i) + 1;
  end
  B(:,i) = bnext;  
  b0(i)  = ym - Xm*b;
  df(i)  = sum(B(:,i)~=0);
  sse(i) = sum((y - X*B(:,i)-b0(i)).^2);
end

stats = struct;
stats.lambda = lambda;
stats.iters = iters;
stats.df = df;
stats.sse = sse;

function [b, active] = cyclicCoordinateDescent(X,y,b,Xss,active,lambda,predictorWeights,XX,Xy)
% Helper function for performing cyclic coordinate descent.
%
active(predictorWeights<=0) = false;
if (any(active))
  % Reindexing the inputs saves some time.
  ind = find(active);
  n = length(active);
  X = X(:,active);
  b = b(active);
  Xss = Xss(active);
  predictorWeights = predictorWeights(active);
  XX = XX(ind,ind);
  Xy = Xy(ind);
  
  r = y - X*b;
  bold = b;
  for i = 1:numel(ind)
    % Covariance update
    % rj'*X(:,i) = <Xj,y> + Xss(j)b(j) - XX_j*b
    b(i) = max(2*(Xy(i) + Xss(i)*b(i)-XX(i,:)*b) - lambda/predictorWeights(i), 0)/(2*Xss(i)); % Soft threshold
    % Update residual to reflect the new value of bj
    r  = r + X(:,i)*(bold(i) - b(i));
  end
  
  bb = b;
  b = 0*(1:n)';
  b(ind) = bb;
  active = (b~=0)';
end

function active = thresholdActive(X,y,b,active,lambda,predictorWeights)
% Returns a list of coefficients to add to the active list.
% 
% The logic is that the potential update to weight j is
% max(2*(rj'*X(:,j)) - lambda/wj,0), where rj is the residual with j
% removed. The idea now is that we're considering the inactive
% predictors, and their partial residual is just the residual
% itself. Hence we're checking to see if any predictors not in the
% model will have non-zero weights after shrinkage.
predictorWeights = abs(predictorWeights); % Should always be positive
active(predictorWeights<=0) = false;
r = y - X(:,active)*b(active);
active = (2*r'*X)>(lambda./predictorWeights');
