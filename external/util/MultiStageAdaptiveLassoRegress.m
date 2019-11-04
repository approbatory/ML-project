function [b,b0,stats,exitFlag] = MultiStageAdaptiveLassoRegress(X,y,varargin)
% function [b,b0,stats,exitFlag] = MultiStageAdaptiveLassoRegress(X,y,varargin)
%
% Uses the adaptive Lasso to prune small weights from the lasso
% solution, by repeatedly applying the adaptive lasso and using the
% previous stage's results as adaptive lasso weights for the next
% stage. 
%
% By default, The number of stages is 2. Setting this to 1 returns the
% lasso solution with smallest cross validation error. Values higher
% than 1 indicate the number of adaptive lasso stages to apply.
%
% See also: LASSOREGRESS, LASSOCROSSVALIDATE

Args = inputParser;
Args.KeepUnmatched = true;
Args.addOptional('numStages',1, @(x) x>0);
Args.addOptional('kfold',8,@(x) x>0);
Args.addOptional('stopOnStable',true);
Args.addOptional('tol',1e-10);
Args.addOptional('LassoFunction',@LassoRegress);
Args.parse(varargin{:});

LASSOF = Args.Results.LassoFunction;
tol = Args.Results.tol;
stopOnStable = Args.Results.stopOnStable;
numStages = Args.Results.numStages;
kfold = Args.Results.kfold;

[numSamples, numRegressors] = size(X);

b = 0*X(1,:)';
b0 = 0;
stats = struct;
stats.B = b*zeros(1,numStages);
stats.b0= zeros(1,numStages);

predictorWeights = ones(numRegressors,1);

exitFlag = 0;
for i = 1:numStages
  [B,b0,lambda] = LASSOF(X,y,'predictorWeights',predictorWeights);
  [indBestLambda,mu,sem] = LassoCrossValidate(X,y,lambda,kfold,0,'predictorWeights',predictorWeights,'LassoFunction',LASSOF);
  b = B(:,indBestLambda);
  b0= b0(indBestLambda);
  
  predictorWeights = b;

  stats.bestLambda(i) = lambda(indBestLambda);
  stats.indBestLambda(i) = indBestLambda;
  stats.lambda(:,i) = lambda;
  stats.mu(:,i) = mu;
  stats.sem(:,i) = sem;
  stats.B(:,i) = b;
  stats.b0(i) = b0;
  stats.df(i) = sum(abs(b)>tol);
  stats.lastStage = i;
  
  if (~any(b))
    exitFlag = 1;
    break;
  elseif (i>1 && (stats.df(i)==stats.df(i-1)) && stopOnStable)
    exitFlag = 2;
    break;
  end
end
  
