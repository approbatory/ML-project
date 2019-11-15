function [frac_correct] = cpd_regress_trial(cpd, meta)

nfactors = length(cpd.lambda);
depvars = {'start','end','turn','strategy','correct','day'};
nvars = length(depvars);
frac_correct = zeros(nfactors,nvars);

for v = 1:length(depvars)
    frac_correct(:,v) = regress_var(cpd, meta, depvars{v});
end

figure()
bar(frac_correct)
legend(depvars)

function [frac_correct] = regress_var(cpd, meta, depvar)
% CPD_PREDICT(cpd, meta, trial_type)
%
% Uses multinomial logistic regression to predict day, strategy, start arm,
% etc. from trial factors.

if strcmp(depvar, 'strategy')
    % remove non-labelled trials
    idx = strcmp('allo-north', meta.strategy) | ...
          strcmp('allo-south', meta.strategy) | ...
          strcmp('ego-right', meta.strategy) | ...
          strcmp('ego-left', meta.strategy);

    X = cpd.factors.trial(idx,:);
    Y = categorical(meta.strategy(idx));
else
    X = cpd.factors.trial;
    Y = categorical(meta.(depvar));
end

nfactors = size(X,2);
frac_correct = zeros(nfactors,1);
for r = 1:nfactors
    B = mnrfit(X(:,r),Y);
    yhat = mnrval(B,X(:,r));
    [~,yi] = max(yhat,[],2);
    yl = categories(Y);
    frac_correct(r) = mean(yl(yi) == Y);
end

% decoding on the full tensor decomposition

% % find leave-one-out error rate
% cv_rate = mean(crossval(@fit_predict,X,Y,'kfold',10));

% % fit the full model, return coefficients
% B = mnrfit(X,Y);

% function frac_correct = fit_predict(xT,yT,xt,yt)
%     % train on (xT,yT) then predict yt using xt
%     B = mnrfit(xT,yT);
%     yhat = mnrval(B,xt);
%     [~,yi] = max(yhat,[],2);
%     yl = categories(yT);
%     frac_correct = mean(yl(yi) == yt);
