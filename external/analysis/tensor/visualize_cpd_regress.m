function [Bn,cv_rate] = visualize_cpd_regress(cpd, meta)
% VISUALIZE_CPD_REGRESS(cpd, meta)
%
% Runs cpd_regress_trial(...) on all dependent variables and
% plots the regression coefficients in a bar chart

% number of factors
nf = length(cpd.lambda);

% dependent variables to visualize
dep_vars = {'start', 'end', 'correct', 'day', 'strategy'};
ndv = length(dep_vars);

% storage for regresion coeffs and cv_rate
B = zeros(nf+1,ndv);
cv_rate = zeros(1,ndv);

for a = 1:ndv
    disp(['Regressing on ',dep_vars{a},'...']);
    [B_, cv_rate(a)] = cpd_regress_trial(cpd, meta, dep_vars{a});
    B(:,a) = sum(abs(B_),2);
end

% normalize regression coefficients
Bn = B ./ repmat( sqrt(sum(B.^2)), size(B,1), 1);

% make plot showing strength of regression coefficients
figure()
bar(0:15,Bn)
set(gca,'Xtick',0:15)
legend({ sprintf('start location (%1.3f accuracy)',cv_rate(1)), ...
         sprintf('end location (%1.3f accuracy)',cv_rate(2)),...
         sprintf('correct choice (%1.3f accuracy)',cv_rate(3)),...
         sprintf('day/session (%1.3f accuracy)',cv_rate(4)),...
         sprintf('internal strategy (%1.3f accuracy)',cv_rate(5))});
xlim([-1 16])
ylabel('abs val regression coefficients (normalized)')
xlabel('trial factors (0 = intercept term)')
